#include "scene.h"

#include "objloader.h"
#include "ray.h"

#include <cmath>
#include <iostream>
#include <cassert>
#include <exception>
#include <random>

#include <optixpp_namespace.h>
#include "sutil.h"


//template std::array<std::vector<ray_physics::segment>, 256> scene::cast_rays<256>();
template std::array<std::array<std::vector<ray_physics::segment>,5>, 512> scene::cast_rays<5,512>(transducer_ & transducer);

scene::scene(const nlohmann::json & config, transducer_ & transducer, optix::Context context) :
    transducer(transducer),
    context(context)
{
    createMaterialPrograms(context, false, closest_hit, any_hit);
    try
    {
        parse_config(config);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error{ "Error while loading scene: " + std::string{ex.what()} };
    }

    //create_empty_world();

    //init();
}

scene::~scene()
{
    destroy_world();
}

/*
//void scene::init()
//{
    for (auto & mesh : meshes)
    {
        //const auto full_path = working_dir + mesh.filename;

        //auto object = add_rigidbody_from_obj(full_path, mesh.deltas, scaling);

        //object->setUserPointer(&mesh);
    }
//}
*/

optix::float3 min_float3(const optix::float3& bbox1, const optix::float3& bbox2 ){
    return optix::make_float3( std::min(bbox1.x, bbox2.x),
                               std::min(bbox1.y, bbox2.y),
                               std::min(bbox1.z, bbox2.z)
                       );
}

optix::float3 min_bbox( std::vector<OptiXMesh> meshes ){
    optix::float3 min_bbox = optix::make_float3(std::numeric_limits<float>::max());
    for (std::vector<OptiXMesh>::iterator it = meshes.begin(); it != meshes.end(); it++){
        min_bbox = min_float3( min_bbox, it->bbox_min );
    }
    return min_bbox;
}

optix::float3 max_float3(const optix::float3& bbox1, const optix::float3& bbox2 ){
    return optix::make_float3(std::max(bbox1.x, bbox2.x),
                              std::max(bbox1.y, bbox2.y),
                              std::max(bbox1.z, bbox2.z)
                       );
}

optix::float3 max_bbox( std::vector<OptiXMesh> meshes ){
    optix::float3 max_bbox = optix::make_float3(std::numeric_limits<float>::min());
    for (std::vector<OptiXMesh>::iterator it = meshes.begin(); it != meshes.end(); it++){
        max_bbox = max_float3( max_bbox, it->bbox_max );
    }
    return max_bbox;
}

template<unsigned int sample_count,unsigned int ray_count>
std::array<std::array<std::vector<ray_physics::segment>,sample_count>, ray_count> scene::cast_rays(transducer_ & transducer)
{
    using namespace ray_physics;

    //arreglo que guarda todos los segmentos de cada sample de cada rayo.
    std::array<std::array<std::vector<segment>,sample_count>, ray_count> segments;

    unsigned int tests = 0;
    unsigned int total_collisions = 0;

    ///step the simulation
    if (m_dynamicsWorld)
    {
        const float ray_start_step { 0.02f };

        for (auto & sample_vector : segments)
        {
            for (auto & segments_vector : sample_vector )
            {
                segments_vector.reserve(ray::max_depth);
            }
        }
        //para cada rayo del trasductor
        //#pragma omp parallel for
        for (size_t ray_i = 0; ray_i < ray_count; ray_i++)
        {
            auto & samples_vector = segments[ray_i];//santi ver si quitar o no
            // modificar xq ahora tengo un arreglo samples. y cada sample tiene un array de segments.

            std::array<ray,sample_count> samples;

            // Add first ray
            {
                ray first_ray
                {
                    //transducer_pos + btVector3(0,0,ray_start_step * ray_i),
                    transducer.element(ray_i).position,                          // from
                            transducer.element(ray_i).direction,                         // initial direction
                            0, //depth                                                          // depth
                            materials.at(starting_material),
                            nullptr,  //material outside
                            initial_intensity / sample_count,//reparto equitativamente la intensidad entre cada uno de los path
                            transducer.frequency,
                            units::length::millimeter_t(0),                           // distance traveled
                            0                                                            // previous ray
                };
                for (size_t sample_i = 0; sample_i < sample_count; sample_i++)
                {
                    samples.at(sample_i) = first_ray;
                }
            }
            for(size_t i = 0; i < ray::max_depth; i++)
            {
                //para cada sample: calculo nuevo segmento y lo guardo. actualizo sample.
                //#pragma omp parallel for
                for (size_t sample_i = 0; sample_i < sample_count; sample_i++)
                {
                    // Pop a ray from the samples array and check if it collides
                    auto & ray_ = samples.at(sample_i);
                    if (!ray_.null)
                    {
                        float r_length = ray_physics::max_ray_length(ray_);
                        auto to = ray_.from + enlarge(ray_.direction, r_length);

                        optix::float3 v = ray_.from + 0.1f * ray_.direction;
                        btCollisionWorld::ClosestRayResultCallback closestResults( btVector3(v.x, v.y, v.z) , btVector3(to.x, to.y, to.z) );

                        m_dynamicsWorld->rayTest(btVector3(v.x, v.y, v.z) , btVector3(to.x, to.y, to.z),closestResults);
                        tests++;

                        if (closestResults.hasHit())
                        {
                            // Substract ray intensity according to distance traveled
                            auto distance_before_hit = ray_.distance_traveled;
                            auto intensity_before_hit = ray_.intensity;

                            const auto organ = static_cast<mesh*>(closestResults.m_collisionObject->getUserPointer());
                            assert(organ);

                            //Simulacion de la penetracion de la onda sobre el tejido para generar interfaces mas difusas y no tan "plasticas"
                            //cada material tiene un parametro thickness (s) que representa la media de los mm que puede penetrar la su superficie la onda.
                            //Por medio de una variable aleatoria calculamos la penetracion de cada rayo sobre la superficie. var: q variable que sigue una probabilidad normal de media 0 y desvio s.
                            std::random_device rd;  //Will be used to obtain a seed for the random number engine
                            std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
                            std::normal_distribution<double> distribution(0.0,organ->material_inside.thickness);
                            float q = std::abs(distribution(generator));

                            //modifico el hit_point para simular la penetracion del rayo - actualizar la intensidad restante que queda al penetrar la malla.

                            auto inside_point = q * ray_.direction + optix::make_float3( closestResults.m_hitPointWorld.getX(), closestResults.m_hitPointWorld.getY(), closestResults.m_hitPointWorld.getZ());
                            ray_physics::travel(ray_, distance_in_mm(ray_.from, inside_point));//aca
                            std::cout << "inside : " << distance_in_mm(ray_.from, inside_point) << std::endl;
                            std::cout << "original" << distance_in_mm(ray_.from, optix::make_float3( closestResults.m_hitPointWorld.getX(), closestResults.m_hitPointWorld.getY(), closestResults.m_hitPointWorld.getZ() ) ) << std::endl;
                            // Calculate refraction and reflection directions and intensities

                            auto result = ray_physics::hit_boundary( ray_, optix::make_float3( closestResults.m_hitPointWorld.getX() , 
									 		       closestResults.m_hitPointWorld.getY() ,
											       closestResults.m_hitPointWorld.getZ()),
                                                                     optix::make_float3( closestResults.m_hitNormalWorld.getX(), 
											 closestResults.m_hitNormalWorld.getY(), 
											 closestResults.m_hitNormalWorld.getZ() ), 
								     *organ);//aca

                            // Register collision creating a segment from the beggining of the ray to the collision point
                            samples_vector[sample_i].emplace_back(segment{ray_.from, inside_point, ray_.direction, result.reflected_intensity, intensity_before_hit, ray_.media.attenuation, distance_before_hit, ray_.media});//aca

                            // Spawn reflection or refraction rays
                            if (result.returned.intensity > ray::intensity_epsilon)
                            {
                                result.returned.parent_collision = samples_vector[sample_i].size()-1;
                                samples.at(sample_i) = result.returned;
                            }
                            else
                                samples.at(sample_i).null = true;

                        }
                        else
                        {
                            // Ray did not reach another media, add a data point at its end.
                            samples_vector[sample_i].emplace_back(segment{ray_.from, to, ray_.direction, 0.0f, ray_.intensity, ray_.media.attenuation, ray_.distance_traveled, ray_.media});
                            samples.at(sample_i).null = true;
                        }
                    }

                }
            }
        }

        for (auto & collision_vector : segments)
        {
            total_collisions += collision_vector.size();
        }
    }

    const float fps = 1.0f/(float( clock() - frame_start ) /  CLOCKS_PER_SEC);
    std::cout << fps << " " << tests << " " << total_collisions << std::endl;
    frame_start = clock();

    return segments;
}

void createMaterialPrograms( optix::Context context, bool use_textures, optix::Program& closest_hit, optix::Program& any_hit)
{
    const char *ptx = sutil::getPtxString( NULL, "phong.cu");
    if ( !closest_hit )
    {
        closest_hit = context->createProgramFromPTXString( ptx,
            use_textures ? "closest_hit_radiance_textured" : "closest_hit_radiance" );
    }
    if ( !any_hit )
        any_hit = context->createProgramFromPTXString( ptx, "any_hit_shadow" );
}

void scene::parse_config(const nlohmann::json & config)
{
    using namespace units::angle;

    working_dir = config.find("workingDirectory") != config.end() ? config.at("workingDirectory") : "";

    const auto & t_pos = config.at("transducerPosition");
    transducer_pos = {t_pos[0], t_pos[1], t_pos[2]};

    const auto & orig = config.at("origin");
    origin = {orig[0], orig[1], orig[2]};

    const auto & spac = config.at("spacing");
    spacing = {spac[0], spac[1], spac[2]};

    starting_material = config.at("startingMaterial");

    scaling = config.at("scaling");

    const auto & mats = config.at("materials");
    if (mats.is_array())
    {
        for (const auto & mat : mats)
        {
            optix::Material m = context->createMaterial();
            m->setClosestHitProgram(0u, closest_hit);
            m->setAnyHitProgram( 1u, any_hit);
            m["impedance"]->setFloat( mat["impedance"].get<float>() );
            m["attenuation"]->setFloat( mat["attenuation"].get<float>() );
            m["mu0"]->setFloat( mat["mu0"].get<float>() );
            m["mu1"]->setFloat( mat["mu1"].get<float>() );
            m["sigma"]->setFloat( mat["sigma"].get<float>() );
            m["specularity"]->setFloat( mat["specularity"].get<float>() );
            m["shininess"]->setFloat( mat["shininess"].get<float>() );
            m["thickness"]->setFloat( mat["thickness"].get<float>() );
            materials[mat["name"].get<std::string>()] = m;
        }
    }
    else
    {
        throw std::runtime_error("materials must be an array");
    }

    const auto & meshes_ = config.at("meshes");
    if (meshes_.is_array())
    {
        optix::GeometryGroup geometry_group = context->createGeometryGroup();
        for (const auto & mesh_i : meshes_)
        {
            const auto & deltas { mesh_i.at("deltas") };
            OptiXMesh m;
            m.context = context;
            m.use_tri_api = true;
            m.ignore_mats = false;
            loadMesh(mesh_i.at("file"), m);
            m.geom_instance["rigid"]->setInt(static_cast<bool>(mesh_i.at("rigid")));
            m.geom_instance["vascular"]->setInt(static_cast<bool>(mesh_i.at("vascular")));
            m.geom_instance["deltas"]->setFloat(optix::make_float3(deltas[0], deltas[1], deltas[2]));
            m.geom_instance["outsideNormals"]->setInt(static_cast<bool>(mesh_i.at("outsideNormals")));  //Outside Material should go in Ray Payload
            std::string outside_material = working_dir + mesh_i["outside_material"].get<std::string>();
            m.geom_instance["outsideMaterial"]->setUserData(sizeof(outside_material), &outside_material);
            geometry_group->addChild(m.geom_instance);
            meshes.emplace_back(m);
        }
        geometry_group->setAcceleration( context->createAcceleration("trbvh") );
        aabb.set(min_bbox(meshes), max_bbox(meshes));
        context["scene"]->set(geometry_group);


    }
    else
    {
        throw std::runtime_error("meshes must be an array");
    }
}

/*void scene::create_empty_world()
{
    m_collisionConfiguration = std::make_unique<btDefaultCollisionConfiguration>();

    m_dispatcher = std::make_unique<btCollisionDispatcher>(m_collisionConfiguration.get());

    m_broadphase = std::make_unique<btDbvtBroadphase>();

    m_solver = std::make_unique<btSequentialImpulseConstraintSolver>();

    m_dynamicsWorld = std::make_unique<btDiscreteDynamicsWorld>(m_dispatcher.get(),m_broadphase.get(),m_solver.get(),m_collisionConfiguration.get());

    m_dynamicsWorld->setGravity(btVector3(0,-10,0));
}*/

void scene::destroy_world()
{
    //delete collision shapes
   /* for (int j = 0; j < m_collisionShapes.size(); j++)
    {
        btCollisionShape* shape = m_collisionShapes[j];
        delete shape;
    }
    m_collisionShapes.clear();

    m_dynamicsWorld.reset();
    m_solver.reset();
    m_broadphase.reset();
    m_dispatcher.reset();
    m_collisionConfiguration.reset();*/

    context->destroy();
}

units::length::millimeter_t scene::distance_in_mm(const optix::float3 & v1, const optix::float3 & v2) const
{
    using namespace std;

    auto x_dist = abs(v1.x - v2.x) * spacing[0];
    auto y_dist = abs(v1.y - v2.y) * spacing[1];
    auto z_dist = abs(v1.z - v2.z) * spacing[2];

    return units::length::millimeter_t(sqrt(pow(x_dist,2) + pow(y_dist,2) + pow(z_dist,2)) * 10);
}

optix::float3 scene::enlarge(const optix::float3 & versor, float mm) const
{
    assert( pow( static_cast<double>( optix::length(versor) ), 2.) < 1.1);
    return mm / 100.0f * optix::make_float3(  spacing[0] * versor.x,
                                            spacing[1] * versor.y,
                                            spacing[2] * versor.z );
}

btRigidBody * scene::add_rigidbody_from_obj(const std::string & fileName, std::array<float, 3> deltas, float scaling)
{
    GLInstanceGraphicsShape* glmesh = load_mesh_from_obj(fileName, "");
    printf("[INFO] Obj loaded: Extracted %d verticed from obj file [%s]\n", glmesh->m_numvertices, fileName.c_str());

    const GLInstanceVertex& v = glmesh->m_vertices->at(0);
    btTriangleIndexVertexArray* tiva = new btTriangleIndexVertexArray(glmesh->m_numIndices / 3, &glmesh->m_indices->at(0), 3* sizeof(int),
                                                                      glmesh->m_numvertices, (btScalar*)(&(v.xyzw[0])), sizeof(GLInstanceVertex));

    btBvhTriangleMeshShape* shape = new btBvhTriangleMeshShape(tiva, true);

    m_collisionShapes.push_back(shape);

    float _scaling[4] = {scaling,scaling,scaling,1};

    btVector3 localScaling(_scaling[0],_scaling[1],_scaling[2]);
    shape->setLocalScaling(localScaling);

    btTransform startTransform;
    startTransform.setIdentity();

    //std::array<float, 3> origin { -18, -22, -5 }; // origin for organs scene
    float pos[4] = {deltas[0]*_scaling[0]*_scaling[0],deltas[1]*_scaling[1]*_scaling[1],deltas[2]*_scaling[2]*_scaling[2],0};
    btVector3 position(pos[0] + origin[0], pos[1] + origin[1], pos[2] + origin[2]);
    startTransform.setOrigin(position);

    btScalar mass(0.f);
    btVector3 localInertia(0, 0, 0);
    auto myMotionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo cInfo(mass, myMotionState, shape, localInertia);
    auto body = new btRigidBody(cInfo);
    body->setUserIndex(-1);
    m_dynamicsWorld->addRigidBody(body);
    return body;
}

void scene::step(float delta_time)
{
    m_dynamicsWorld->stepSimulation(delta_time);
}

// TODO: Is this equals to distance_in_mm?
units::length::millimeter_t scene::distance(const optix::float3 & from, const optix::float3 & to) const
{
    // TODO: Use scaling in this calculation
    optix::float3 v = to - from;
    return units::length::millimeter_t( optix::length(v) * 10.0f );
}


