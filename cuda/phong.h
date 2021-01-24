/* 
 * Copyright (c) 2018, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
//#include <cmath>
#include <math_functions.h>
#include <optix_world.h>
#include "common.h"
#include "helpers.h"
#include <optix_device.h>
#include "random.h"

using namespace optix;

/*
struct PerRayData_radiance
{
    float3 result;
    float importance;
    int depth;
    float media_impedance;
    float media_attenuation;
};

struct PerRayData_shadow
{
    float3 attenuation;
};*/


struct PRD{
    float3 origin; //hit_point antes. Origen del rayo. Podria no usar mas rtCurrentRay
    float3 direction; //Direccion del rayo
    float importance; //intensity
    unsigned int depth;
    float distance_traveled;//units::length::millimeter_t distance_traveled;
    float time_traveled;
    float media_impedance;
    float media_attenuation;
    float media_mu0;
    float media_mu1;
    float media_sigma;
    float btti; //Back-To-Transducer Intensity
    int bounce_or_transmission; //0=bounce, 1=transmission
    int sample;
    int mesh_index;
    float is_vascular;
};


rtDeclareVariable(uint,              max_depth, , );
rtDeclareVariable(float3,            ambient_light_color, , );
rtDeclareVariable(float,             scene_epsilon, , );
rtDeclareVariable(rtObject,          top_object, , );
rtDeclareVariable(rtObject,          top_shadower, , );
rtDeclareVariable(uint,              max_rows,,);
rtDeclareVariable(float3,            spacing,,);
rtDeclareVariable(float,             transducer_frequency,,);
rtDeclareVariable(float,             speed_of_sound,,);
rtDeclareVariable(unsigned int,      samples_te,,);
rtDeclareVariable(unsigned int,      transducer_elements,,);
rtBuffer<float, 3>                   results;
rtBuffer<float>                      normal_dist;
rtBuffer<float>                      uniform_dist;


//rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float, t_hit, rtIntersectionDistance, );
rtDeclareVariable(PRD, prd, rtPayload, );
rtDeclareVariable(uint, launch_id, rtLaunchIndex,);

//rtDeclareVariable(PerRayData_shadow,   prd_shadow, rtPayload, );

static __device__ void phongShadowed()
{
    // this material is opaque, so it fully attenuates all shadow rays
    //prd_shadow.attenuation = optix::make_float3(0.0f);
    //rtTerminateRay();

    rtIgnoreIntersection();
    rtTerminateRay();//rtIgnoreIntersection();
}

static __device__ float power_cosine_variate(int v){
    //std::random_device rd;
    //std::mt19937 generator(rd());
    //std::uniform_real_distribution<double> distribution(0.0,1.0);
    //double number = distribution(generator);
    unsigned int index = launch_id + transducer_elements * (prd.depth + max_depth * prd.sample) + 1;
    index = index % (max_depth*transducer_elements*samples_te);
    float number = uniform_dist[index];
    int idx = v + 1;
    float exp = (double)1.0 / idx;
    return pow(number,exp);
}

static __device__ float3 random_unit_vector(optix::float3 v, float cos_theta)
{
    bool flag = false;
    float px, py,p;
    unsigned int next_try = 0;
    unsigned int index = 0;
    do
    {

        //Generating a random point within a circle (uniformly)
        //https://programming.guide/random-point-within-circle.html
        //std::random_device rd;  //Will be used to obtain a seed for the random number engine
        //std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
        //std::uniform_real_distribution<double> distribution(0.0,1.0);
        index = (launch_id + transducer_elements * (prd.depth + max_depth * prd.sample) + next_try) % (max_depth*samples_te*transducer_elements);

        float rand1 = uniform_dist[index];
        float rand2 = uniform_dist[index*2 % (max_depth*samples_te*transducer_elements)];
        float a = rand1 * 2 * M_PIf;
        float r = 0.5 * sqrtf(rand2);

        //in Cartesian coordinates
        px = r * cos((float)a);
        py = r * sin((float)a);
        p = px * px + py * py;
        next_try++;
    } while (! (p <= 0.25) );
    float vx = v.x;
    float vy = v.y;
    float vz = v.z;
    if ( abs(vx) > abs(vy) )
    {
        vx = vy;
        vy = v.x;
        flag = true;
    }
    float b = 1 - vx * vx;
    float radicando = 1 - cos_theta * cos_theta;
    radicando = radicando / (p * b);
    float c = sqrt(radicando);
    px = px * c;
    py = py * c;
    float d = cos_theta - vx * px;
    float wx = vx * cos_theta - b * px;
    float wy = vy * d + vz * py;
    float wz = vz * d - vy * py;
    if (flag)
    {
        float aux = wy;
        wy = wx;
        wx = aux;
    }
    return optix::make_float3(wx,wy,wz);
}

static __device__ float3 snells_law(float3 direction, float3 normal, float incidence_angle, float refraction_angle, float refraction_ratio){

    const optix::float3 & l = direction;
    const optix::float3 & n = normal;
    const float c = incidence_angle;
    const float r = refraction_ratio;

    return (r * l + (r*c - sqrtf(1-r*r*(1-c*c))) * n );

}

static __device__ float reflection_intensity(float intensity_in, float media_1, float incidence_angle, float media_2, float refracted_angle)
{

    float num = media_1 * incidence_angle - media_2 * refracted_angle;
    float denom = media_1 * incidence_angle + media_2 * refracted_angle;
    float factor = ( num / denom ) * ( num / denom );

    return intensity_in * factor;
}

static __device__ float reflected_intensity(const optix::float3 direction, const optix::float3 refraction_direction, const optix::float3 reflection_direction, const float & specularity)
{

    // Eq. 8 in Mattausch
    float dist = sqrt(dot(direction,direction));

   float refraction_angle = optix::dot( direction, refraction_direction ) / dist; //cos theta1?
   const float refraction_factor = powf(refraction_angle, specularity);
   float reflection_angle = optix::dot( direction, reflection_direction ) / dist; //cos theta2?
   const float reflection_factor = powf(reflection_angle, specularity);
   //std::cout << "refraction factor : " << refraction_factor << " reflection factor : " << reflection_factor << std::endl;
   return fmaxf(refraction_factor,0.0f) + fmaxf(reflection_factor,0.0f);

}

float distance(optix::float3 & v1, optix::float3 & v2)
{

    float x_dist = abs(v2.x - v1.x) * spacing.x;
    float y_dist = abs(v2.y - v1.y) * spacing.y;
    float z_dist = abs(v2.z - v1.z) * spacing.z;

    return sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist) * 10;
}

void travel( float distance, float attenuation ){

    prd.distance_traveled += distance;
    prd.importance = prd.importance * expf(-attenuation * distance * 0.01 * transducer_frequency );
                                                                      //That 0.01 should be 0.1
}

// impedance, attenuation, mu0, mu1, sigma, specularity, shininess, thickness
static
__device__ void phongShade( float  impedance,
                            float  attenuation,
                            float  mu0,
                            float  mu1,
                            float  sigma,
                            float  specularity,
                            float  shininess,
                            float  thickness,
                            float3 ffnormal,
                            int index,
                            int vascular)
{

    //OBTENER NUMERO RANDOM PARA SIMULAR PENETRACION DE RAYO
    //Generador de numeros aleatorios --> Distribucion normal --> Media=0.0, Desv.Est.=thickness
    //Modificar hit_point de acuerdo al numero random

    prd.time_traveled = t_hit;

    float x = normal_dist[launch_id + transducer_elements * (prd.depth + max_depth * prd.sample)]*thickness;
    //rtPrintf("X: %f\n",x);
    float3 hit_point = prd.origin + t_hit * prd.direction;
    float3 new_hit_point = hit_point + x * prd.direction;
    float3 v = new_hit_point - prd.origin;  //inside_point


    travel( distance(prd.origin,   new_hit_point), prd.media_attenuation );
    //prd.distance_traveled = distance(prd.origin, hit_point);

    //Calcular distancia recorrida.
    //float dist = sqrt(dot(v,v));

    float random_angle = power_cosine_variate(shininess);

    float3 random_normal = random_unit_vector(ffnormal, random_angle);
    //rtPrintf("Normal: %f,%f,%f\nRandom normal: %f,%f,%f - Angle: %f\n", ffnormal.x,ffnormal.y,ffnormal.z,random_normal.x,random_normal.y,random_normal.z,random_angle);
    float incidence_angle = optix::dot( v, -random_normal); // cos theta_1
    if (incidence_angle < 0)
        incidence_angle = optix::dot( v, random_normal);

    const float refr_ratio = prd.media_impedance / impedance;

    float refraction_angle = 1 - refr_ratio*refr_ratio * (1 - incidence_angle*incidence_angle);
    const bool total_internal_reflection = refraction_angle < 0;

    //if ((impedance == 7.8f) && (prd.media_impedance == 1.65f))
    //    rtPrintf("incidence_angle=%f - refr_ratio=%f\n", incidence_angle, refr_ratio);
    refraction_angle = sqrt(refraction_angle);

    float3 refraction_direction = snells_law(prd.direction, random_normal, incidence_angle, refraction_angle, refr_ratio);
    refraction_direction = normalize(refraction_direction);

    float3 reflection_direction = prd.direction + 2*incidence_angle*random_normal;
    reflection_direction = normalize(reflection_direction);



    float intensity_reflected = 0.f;
    if (total_internal_reflection) {
        intensity_reflected = prd.importance;
    }
    else{
        intensity_reflected = reflection_intensity(prd.importance,prd.media_impedance,incidence_angle,impedance,refraction_angle);
    }


    const float intensity_refracted = prd.importance - intensity_reflected;

    const float back_to_transducer_intensity = reflected_intensity(v, refraction_direction, reflection_direction, specularity) * random_angle;
    float reflection_probability = intensity_reflected / prd.importance;
    prd.origin = new_hit_point;

    x = uniform_dist[launch_id + transducer_elements * (prd.depth + max_depth * prd.sample)];
    //if (launch_id==0){rtPrintf("Launch_id: %d - transducer_elements: %d - depth: %d - max_depth: %d - sample: %d - max_samples: %d\n", launch_id, transducer_elements, prd.depth, max_depth, prd.sample, samples_te);
    //rtPrintf("X: %d, %f\n",launch_id + transducer_elements * (prd.depth + max_depth * prd.sample),x);}

    prd.btti = back_to_transducer_intensity;

    if ( ( reflection_probability > x ) )
    {

        prd.direction = reflection_direction;
        prd.importance = intensity_reflected > scene_epsilon ? intensity_reflected : 0.0f;
        prd.bounce_or_transmission = 0;
        prd.is_vascular = vascular != 0 ? 1.f : 0.f;
        //rtPrintf("is_vascular: %f\n", prd.is_vascular);


    }
    else {

        prd.direction = refraction_direction;
        prd.importance = intensity_refracted > scene_epsilon ? intensity_refracted : 0.0f;
        prd.bounce_or_transmission = 1;
        prd.mesh_index = index;
        prd.media_attenuation = attenuation;
        prd.media_impedance = impedance;
        prd.media_mu0 = mu0;
        prd.media_mu1 = mu1;
        prd.media_sigma = sigma;
        prd.is_vascular = vascular != 0 ? 1.f : 0.f;
        //rtPrintf("is_vascular: %f\n", prd.is_vascular);
    }



}

