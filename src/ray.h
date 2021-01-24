#ifndef RAY_H
#define RAY_H

//#include <LinearMath/btVector3.h>

#include <units.h>
#include "mesh.h"
#include <optixu_math_namespace.h>
//class material;
class mesh;

namespace ray_physics {

struct ray          //Reemplazar con rayo optix
{
    optix::float3 from, direction;
    size_t depth;
    material media;
    const material * media_outside = nullptr; // keep track of media outside of vascularities, so we can switch back to it when the ray gets out of them
    float intensity, frequency;
    units::length::millimeter_t distance_traveled; // [mm]
    unsigned short parent_collision; // position in collision vector

    static constexpr size_t max_depth = 10;
    static constexpr float intensity_epsilon = 1e-10;
    bool null = false;
};

struct segment
{
    optix::float3 from, to, direction;
    float reflected_intensity; // reflected back to the transducer, at the end of the segment
    float initial_intensity, attenuation, mu0, mu1,sigma;
    float distance_traveled; //[mm]
    unsigned int ray, depth;
    int sample;
    bool null = true;
    //units::length::millimeter_t distance_traveled; // traveled from the transducer to the beginning of the segment
    //const material & media;
};

struct collision
{
    optix::float3 position;
    unsigned short parent_collision; // position in collision vector
};

struct hit_result { float reflected_intensity; ray returned; };

hit_result hit_boundary(const ray & r, const optix::float3 &hit_point, const optix::float3 & surface_normal, const mesh & collided_mesh);

// Advance through homogeneous media and decrease intensity accordingly
void travel(ray & r, units::length::millimeter_t mm);

bool should_travel(const ray & r);

float max_ray_length(const ray & r);

optix::float3 snells_law(const optix::float3 & ray_direction, const optix::float3 & surface_normal, float incidence_angle, float refraction_angle, float refr_ratio);

/**
 * Intensity of the reflected ray.
 * IMPORTANT: This it NOT the intensity reflected back to the transducer.
 */
float reflection_intensity(const float intensity_in, const float media_1, const float incidence_angle, const float media_2, const float refracted_angle);

// Intensity reflected back to the transducer when a ray passes through an interface.
float reflected_intensity(const float ray_intensity, const float incidence_angle, const material & ray_media, const material & colliding_media);

float reflected_intensity(const optix::float3 direction, const optix::float3 refraction_direction, const optix::float3 reflection_direction, const material & colliding_media);

optix::float3 random_unit_vector(optix::float3 v, float cos_theta);

float power_cosine_variate(int v);

} // end ray_physics

#endif // RAY_H
