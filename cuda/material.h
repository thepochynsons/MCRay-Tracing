#include <optix_world.h>
#include "../src/mesh.h"
//Luces

struct PerRayData_radiance {
    int depth;
    float intensity;
    float back_to_transducer_intensity;
    int current_material_index;
    int outside_material_index;
    //float transducer_freq; //Dejarlo aca para potencia de 2? o cambiar a variable declarada?
};

rtBuffer<TransducerElement>             transducer_elements; //luces, generadores de rayos.

rtDeclareVariable( int,                 max_depth,          ,                       );
rtDeclareVariable( float,               scene_epsilon,      ,                       );
rtDeclareVariable( float,               intensity_epsilon,  ,                       );

//es necesario?
rtDeclareVariable( int,                 sample_count,       ,                       );

rtDeclareVariable( rtObject,            top_object,         ,                       );
rtDeclareVariable( optix::Ray,          ray,                rtCurrentRay,           );
rtDeclareVariable( float,               t_hit,              rtIntersectionDistance, );
rtDeclareVariable( PerRayData_radiance, prd,                rtPayload,              );

static __device__ void compute_ray(  )
{
    float3 hit_point = ray.origin + t_hit * ray.direction;

    //actualizacion de hit_point para simular penetracion del rayo.
    std::random_device rd;
    std::mt19937 generator( rd( ) );
    //thickness del material interno de la mesh.
    std::normal_distribution<double> distribution( 0.0, thickness );
    float q = std::abs( distribution( generator ) );

    hit_point = hit_point + q;

    const material * material_after_vascularities = nullptr;
    const auto & material_after_collision = [&material_after_vascularities]() -> const material &
    {

    }


}
