#ifndef SCENE_H
#define SCENE_H

//#include "btBulletDynamicsCommon.h"

#include "ray.h"
#include "mesh.h"
#include "transducer.h"

#include <ctime>
#include <memory>
#include <unordered_map>
#include <array>
#include <vector>

#include <nlohmann/json.hpp>
#include <units/units.h>

#include "OptiXMesh.h"


class scene
{
    using transducer_ = transducer<512>;

public:
    explicit scene(const nlohmann::json & config, transducer_ & transducer, optix::Context context);
    ~scene();

    void init();

    template<unsigned int sample_count,unsigned int ray_count>
    std::array<std::array<std::vector<ray_physics::segment>,sample_count>, ray_count>cast_rays(transducer_ & transducer);

    void step(float delta_time);

    units::length::millimeter_t distance(const optix::float3 & from, const optix::float3 & to) const;
    //void setTransducer (transducer_  transducer);

protected:
    std::string working_dir;

    std::unordered_map<std::string, optix::Material> materials;
    std::string starting_material;

    std::vector<OptiXMesh> meshes;

    transducer_ & transducer;

    const float intensity_epsilon { 1e-8 };
    const float initial_intensity { 1.0f };

    std::array<float,3> spacing;
    std::array<float,3> origin;
    float scaling;

    void parse_config(const nlohmann::json & config);
    void createMaterialPrograms( optix::Context context, bool use_textures, optix::Program& closest_hit, optix::Program& any_hit);
    void create_empty_world();
    void destroy_world();

    optix::float3 min_bbox(std::vector<OptiXMesh>);
    optix::float3 max_bbox(std::vector<OptiXMesh>);

    units::length::millimeter_t distance_in_mm(const optix::float3 & v1, const optix::float3 & v2) const;
    optix::float3 enlarge(const optix::float3 & versor, float mm) const;

    class btRigidBody * add_rigidbody_from_obj(const std::string & fileName, std::array<float, 3> deltas, float scaling);

    btAlignedObjectArray<btCollisionShape*>	m_collisionShapes;
    std::unique_ptr<btBroadphaseInterface> m_broadphase;
    std::unique_ptr<btCollisionDispatcher> m_dispatcher;
    std::unique_ptr<btConstraintSolver> m_solver;
    std::unique_ptr<btDefaultCollisionConfiguration> m_collisionConfiguration;
    std::unique_ptr<btDiscreteDynamicsWorld> m_dynamicsWorld;

    btVector3 transducer_pos;
    std::array<units::angle::degree_t, 3> transducer_dir;

    optix::Context context;
    optix::Program closest_hit;
    optix::Program any_hit;
    optix::Aabb    aabb;

    clock_t frame_start;
};

#endif // SCENE_H
