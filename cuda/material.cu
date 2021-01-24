#include <optix.h>
#include <optixu/optixu_math_namespace.h>

using namespace optix;

/*rtDeclareVariable(float, impedance, , );
rtDeclareVariable(float, attenuation_factor, , );
rtDeclareVariable(float, mu0, , );
rtDeclareVariable(float, mu1, , );
rtDeclareVariable(float, sigma, , );
rtDeclareVariable(float, specularity, , );
rtDeclareVariable(float, shininess, , );
rtDeclareVariable(float, thickness, , );*/

rtDeclareVariable(float3, geometric_normal, attribute geometric_normal, );
rtDeclareVariable(float3, shading_normal, attribute shading_normal, );

RT_PROGRAM void closest_hit()
{
    float3 world_shading_normal = normalize( rtTransformNormal( RT_OBJECT_TO_WORLD, shading_normal) );
    float3 world_geometric_normal = normalize( rtTransformNormal( RT_OBJECT_TO_WORLD, geometric_normal) );

    float3 ffnormal = faceforward( world_shading_normal, -ray.direction, world_geometric_normal );
    compute_ray( );


}
