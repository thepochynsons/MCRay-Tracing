#ifndef CAMERA_H
#define CAMERA_H

//#include "btBulletDynamicsCommon.h"

#include "ray.h"
#include "mesh.h"
#include "transducer.h"

#include <ctime>
#include <memory>
#include <unordered_map>
#include <array>
#include <vector>

#include <json.hpp>
#include <units.h>
#include <sutil/sutil.h>
#include "OptiXMesh.h"

#include <optixpp_namespace.h>
//#include <optix_device.h>

//#include <cuda.h>
//#include <curand_kernel.h>
//#include <curand.h>
//#include "device_launch_parameters.h"

using transducer_ = transducer<512>;

/*-------Bullet Vars-------*/
/*
//btAlignedObjectArray<btCollisionShape*>             m_collisionShapes;
//std::unique_ptr<btBroadphaseInterface>              m_broadphase;
//std::unique_ptr<btCollisionDispatcher>              m_dispatcher;
//std::unique_ptr<btConstraintSolver>                 m_solver;
//std::unique_ptr<btDefaultCollisionConfiguration>    m_collisionConfiguration;
//std::unique_ptr<btDiscreteDynamicsWorld>            m_dynamicsWorld;
*/

/*-----Transducer Vars-----*/
optix::float3                           transducer_pos;
std::array<units::angle::degree_t, 3>   transducer_dir;

/*-------OptiX  Vars-------*/
//optix::Context context;
optix::Program closest_hit;
optix::Program any_hit;
optix::Aabb    aabb;

clock_t frame_start;

/*-------Scene Vars-------*/
std::string         working_dir;
optix::float3       spacing;
optix::float3       origin;
float               scaling;
const float         intensity_epsilon { 1e-8f };
const float         initial_intensity { 1.0f };
std::string         starting_material;

std::unordered_map<std::string, optix::Material> materials;
std::vector<OptiXMesh> meshes;

//transducer_ & transducer;

void createMaterialPrograms(
        optix::Context  context,
        bool            use_textures,
        optix::Program& closest_hit,
        optix::Program& any_hit)
{
    const char *ptx = sutil::getPtxString( NULL, "phong.cu");
    if ( !closest_hit )
    {
        closest_hit = context->createProgramFromPTXString( ptx, "closest_hit_radiance" );
    }
    try {
        closest_hit->validate();
    } catch (optix::Exception ex) {
        std::cout << ex.what() << std::endl;
    }

}

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


void parse_config(const nlohmann::json & config, optix::Context context)
{
    using namespace units::angle;

    working_dir = config.find("workingDirectory") != config.end() ? config.at("workingDirectory") : "";

    const auto & t_pos = config.at("transducerPosition");
    transducer_pos = optix::make_float3(t_pos[0], t_pos[1], t_pos[2]);//optix::make_float3(-18., -22., -5.);//

    const auto & orig = config.at("origin");
    origin = optix::make_float3(orig[0], orig[1], orig[2]);
    context["origin"]->setFloat(origin);

    const auto & spac = config.at("spacing");
    spacing = optix::make_float3(spac[0], spac[1], spac[2]);
    context["spacing"]->setFloat(spacing);

    starting_material = config.at("startingMaterial");

    scaling = config.at("scaling");
    int mesh_id = 1;
    const auto & mats = config.at("materials");
    if (mats.is_array())
    {
        for (const auto & mat : mats)
        {
            optix::Material m = context->createMaterial();
            m->setClosestHitProgram(0u, closest_hit);
            m["impedance"]->setFloat( mat["impedance"].get<float>() );
            m["attenuation"]->setFloat( mat["attenuation"].get<float>() );
            m["mu0"]->setFloat( mat["mu0"].get<float>() );
            m["mu1"]->setFloat( mat["mu1"].get<float>() );
            m["sigma"]->setFloat( mat["sigma"].get<float>() );
            m["specularity"]->setFloat( mat["specularity"].get<float>() );
            m["shininess"]->setFloat( mat["shininess"].get<float>() );
            m["thickness"]->setFloat( mat["thickness"].get<float>() );
            materials[mat["name"].get<std::string>()] = m;
            try { m->validate(); }
            catch (optix::Exception ex){
                std::cout << ex.getErrorString() << std::endl;
            }
            if (mat["name"].get<std::string>() == starting_material){
                context["start_attenuation"]->setFloat((float)mat["attenuation"]);
                context["start_impedance"]->setFloat((float)mat["impedance"]);
                context["start_mu0"]->setFloat((float)mat["mu0"]);
                context["start_mu1"]->setFloat((float)mat["mu1"]);
                context["start_sigma"]->setFloat((float)mat["sigma"]);
                //context["start_index"]->setInt((int)m["index"]);
                //printf("Values:\n\t Attenuation: %f\n\t Impedance: %f\n\t Mu0: %f\n\t Mu1: %f\n\t Sigma: %f\n",mat["attenuation"].get<float>(),mat["impedance"].get<float>(),mat["mu0"].get<float>(),mat["mu1"].get<float>(),mat["sigma"].get<float>());
                //std::cout << "Start impedance: " << m["impedance"]->getFloat() << std::endl << "Start attenuation: " << m["attenuation"]->getFloat() << std::endl;
            }
        }
    }
    else
    {
        throw std::runtime_error("Materials must be an array");
    }

    std::cout << "Initial impedance: " << context["start_impedance"]->getFloat() << std::endl << "Initial attenuation: " << context["start_attenuation"]->getFloat() << std::endl;
    std::cout << "Initial mu0: " << context["start_mu0"]->getFloat() << std::endl << "Initial mu1: " << context["start_mu1"]->getFloat() << std::endl << "Initial sigma: " << context["start_sigma"]->getFloat() << std::endl;
    const auto & meshes_ = config.at("meshes");
    if (meshes_.is_array())
    {
        optix::GeometryGroup geometry_group = context->createGeometryGroup();
        for (const auto & mesh_i : meshes_)
        {
            const auto & deltas { mesh_i.at("deltas") };
            //std::cout << "About to load OptiX Mesh" << std::endl;
            OptiXMesh m;
            m.context = context;
            m.use_tri_api = true;
            m.ignore_mats = false;
            m.material = materials[mesh_i.at("material")];

            float x = static_cast<float>(deltas[0])*scaling*scaling + origin.x;
            float y = static_cast<float>(deltas[1])*scaling*scaling + origin.y;
            float z = static_cast<float>(deltas[2])*scaling*scaling + origin.z;

            optix::Matrix4x4 traslation = optix::Matrix4x4::translate(optix::make_float3(x,y,z));
            optix::Matrix4x4 scale = optix::Matrix4x4::scale(optix::make_float3(scaling));
            optix::Matrix4x4 accum = traslation * scale;
            //if ((mesh_i.at("file") == "liver.obj") || (mesh_i.at("file") == "skin.obj") || (mesh_i.at("file") == "right_kidney.obj"))
            //    std::cout << mesh_i.at("file") << std::endl << x<<","<<y<<","<<z << std::endl;
            loadMesh(working_dir + mesh_i["file"].get<std::string>(), m, accum);
            std::cout << mesh_i["file"].get<std::string>() << " - " << mesh_id << std::endl;
            m.geom_instance["index"]->setInt(mesh_id++);
            m.geom_instance["rigid"]->setInt(static_cast<bool>(mesh_i.at("rigid")));
            m.geom_instance["vascular"]->setInt(static_cast<bool>(mesh_i.at("vascular")));
            m.geom_instance["deltas"]->setFloat(optix::make_float3(deltas[0], deltas[1], deltas[2]));
            m.geom_instance["outsideNormals"]->setInt(static_cast<bool>(mesh_i.at("outsideNormals")));

//            std::string outside_material = mesh_i["outsideMaterial"].get<std::string>();
            m.geom_instance->setMaterial(0u, materials[mesh_i["material"]]);
            meshes.emplace_back(m);
            geometry_group->addChild(m.geom_instance);

            try {
                m.geom_instance->validate();
            } catch (optix::Exception e) {
                std::cout << e.getErrorString() << std::endl;
            }

        }
        //aabb.set(min_bbox(meshes), max_bbox(meshes));
        geometry_group->setAcceleration( context->createAcceleration("trbvh") );

        context["top_object"]->set(geometry_group);
        //context["top_shadower"]->set(geometry_group);
        try {
            geometry_group->validate();
            geometry_group->getAcceleration()->validate();
        } catch (optix::Exception e) {
            std::cout << e.getErrorString() << std::endl;
        }



    }
    else
    {
        throw std::runtime_error("Meshes must be an array");
    }
}


void loadScene(const nlohmann::json & config, optix::Context ctx){
    createMaterialPrograms(ctx, false, closest_hit, any_hit);
    ctx["scene_epsilon"]->setFloat(intensity_epsilon);
    try {
        parse_config(config, ctx);
    } catch (const std::exception & ex) {
        throw std::runtime_error{ "Error while loading scene: " + std::string{ ex.what() } };
    }
    //std::cout << "Scene loaded" << std::endl;
}

void destroyScene(optix::Context context){
    context->destroy();
};



units::length::millimeter_t distance(const optix::float3 & from, const optix::float3 & to) {
    // TODO: Use scaling in this calculation
    return units::length::millimeter_t( optix::length(to - from) *10.f) ;
}

//void setTransducer (transducer_  transducer);










//units::length::millimeter_t distance_in_mm(const optix::float3 & v1, const optix::float3 & v2) const;
//optix::float3 enlarge(const optix::float3 & versor, float mm) const;

//class btRigidBody * add_rigidbody_from_obj(const std::string & fileName, std::array<float, 3> deltas, float scaling);



#endif // SCENE_H
