#include <optix_world.h>
#include <optix_device.h>
//#include <units.h>
//#include <curand_uniform.h>
#include "transducer.h"
#include "helpers.h"
#include "ray.h"

//#define MAX_STACK_SIZE 10
#define EMPTY_STACK -1

using namespace optix;

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



using namespace units::literals;

rtDeclareVariable(uint,     max_depth, , "Max amount of bounces");
rtDeclareVariable(uint,     max_rows, , "Max number of rows in output_buffer");
rtDeclareVariable(float,    scene_epsilon, , "Minimal importance admited.");
rtDeclareVariable(float,    axial_resolution, , );
rtDeclareVariable(float,    speed_of_sound, , );
rtDeclareVariable(float3,   bad_color, , "Exception color.");
rtDeclareVariable(rtObject, top_object, ,"Scene.");
rtDeclareVariable(float,    transducer_frequency, , "Transducer frequency.");
rtDeclareVariable(float,    start_attenuation,,"Initial material's attenuation.");
rtDeclareVariable(float,    start_impedance,,"Initial material's impedance.");
rtDeclareVariable(float,    start_mu0,, "Initial material's mu0.");
rtDeclareVariable(float,    start_mu1,, "Initial material's mu1.");
rtDeclareVariable(float,    start_sigma,, "Initial material's sigma.");
rtDeclareVariable(int,      start_index,, "Initial material's index.");
rtDeclareVariable(uint,     launch_id, rtLaunchIndex,);
rtDeclareVariable(uint,     samples_te,, "How many samples per ray are casted.");
rtDeclareVariable(float,    max_travel_time,, );
rtDeclareVariable(uint,     transducer_elements,,);
rtDeclareVariable(float3,   spacing,,);
rtBuffer<float3> transducer_positions_buffer;   //Buffer with Transducer Elements positions
rtBuffer<float3> transducer_directions_buffer;  //Buffer with Transducer Elements directions
rtBuffer<int> results;  //Buffer with distances and BTTIs


rtBuffer<float4, 2> output_buffer;
rtBuffer<ray_physics::segment> segments_buffer;



RT_PROGRAM void raygen(){
    float3 te_pos = transducer_positions_buffer[launch_id];
    float3 te_dir = transducer_directions_buffer[launch_id];
    //if (launch_id == 255u)
    //    rtPrintf("Pos: %f,%f,%f - Dir: %f, %f, %f\n", te_pos.x,te_pos.y,te_pos.z,te_dir.x,te_dir.y,te_dir.z);
    unsigned int sample = 0;
    unsigned int index;
    float r_length;
    float3 to;
    do{
        PRD prd;
        prd.origin = te_pos;
        prd.direction = te_dir;
        prd.importance = 1.f; //attenuation
        prd.depth = 0;
        prd.sample = sample;
        prd.distance_traveled = 0; //[mm]
        prd.time_traveled = 0;
        prd.media_attenuation = start_attenuation; //
        prd.media_impedance = start_impedance; //
        prd.media_mu0 = start_mu0;
        prd.media_mu1 = start_mu1;
        prd.media_sigma = start_sigma;
        prd.btti = 0.f;
        prd.mesh_index = EMPTY_STACK;

        float4 attenuation_impedance_stack[10]; //.x -> attenuation | .y -> impedance | .z -> is_vascular | .w -> material_index
        float4 mu_sigma_stack[10];              //.x -> mu0 | .y -> mu1 | .z -> sigma | .w -> material_index
        int stack_index = 0;
        attenuation_impedance_stack[stack_index] = make_float4(start_attenuation,start_impedance,0.f,EMPTY_STACK);


        PRD old_prd = prd;
        for(;;)         //Iteratively trace rays. When the loop ends, its next bounce will be traced.
        {
            //Check constraints.
            if ((prd.depth >= max_depth) || (prd.importance < scene_epsilon ))
                break;
            //Save old prd to posterior ray marching (Check if assign is Copy or Reference)
            old_prd = prd;
            //Create ray and trace it
            Ray ray = make_Ray(prd.origin, prd.direction, 0, scene_epsilon, RT_DEFAULT_MAX);
            r_length = 10.f /*<- cm to mm*/ * log(scene_epsilon/prd.importance) / attenuation_impedance_stack[stack_index].x * transducer_frequency;
            to = prd.origin + r_length / 100.f * make_float3(spacing.x * prd.direction.x,
                                                             spacing.y * prd.direction.y,
                                                             spacing.z * prd.direction.z); //Para guardar el segmento si no choca con nada mas


            rtTrace(top_object, ray, prd);
            results[launch_id] = results[launch_id]+1;



            index = launch_id + transducer_elements * (prd.depth + max_depth * sample);
            //rtPrintf("Index: %d\n", index);


            if ((old_prd.origin.x == prd.origin.x) && (old_prd.origin.y == prd.origin.y) && (old_prd.origin.z == prd.origin.z)){  //chequea si no se actualizó el origen, quiere decir que no se chocó con nada. Quizas podria estar en un programa MISS
                segments_buffer[index] = ray_physics::segment{       //creo el segmento
                        old_prd.origin,
                        to,//old_prd.origin + max_travel_time * old_prd.direction,
                        old_prd.direction,
                        0.0,//prd.btti,
                        0.0,//old_prd.importance,
                        attenuation_impedance_stack[stack_index].x,
                        mu_sigma_stack[stack_index].x,
                        mu_sigma_stack[stack_index].y,
                        mu_sigma_stack[stack_index].z,
                        old_prd.distance_traveled,
                        launch_id,
                        prd.depth,
                        prd.sample,
                        false
            };
            } else {
                //if (prd.is_vascular == 0.f){
                segments_buffer[index] = ray_physics::segment{       //creo el segmento
                        old_prd.origin,
                        prd.origin,
                        old_prd.direction,
                        prd.btti,
                        old_prd.importance,
                        attenuation_impedance_stack[stack_index].x,
                        mu_sigma_stack[stack_index].x,
                        mu_sigma_stack[stack_index].y,
                        mu_sigma_stack[stack_index].z,
                        old_prd.distance_traveled,
                        launch_id,
                        prd.depth,
                        prd.sample,
                        false
                };/*} else {
                    segments_buffer[index] = ray_physics::segment{       //creo el segmento
                            old_prd.origin,
                            prd.origin,
                            old_prd.direction,
                            0.1f, //para definir manualmente el brillo de las paredes de los vasos
                            old_prd.importance,
                            attenuation_impedance_stack[stack_index].x,
                            mu_sigma_stack[stack_index].x,
                            mu_sigma_stack[stack_index].y,
                            mu_sigma_stack[stack_index].z,
                            old_prd.distance_traveled,
                            launch_id,
                            prd.depth,
                            prd.sample,
                            false
                    };
                }*/
            }



            //Update depth value

            if ((prd.bounce_or_transmission) && (stack_index >= 0)){
                if ((prd.mesh_index == attenuation_impedance_stack[stack_index].w) && (prd.mesh_index != EMPTY_STACK)){ //choco el mismo mesh, entonces estoy saliendo
                    attenuation_impedance_stack[stack_index] = make_float4(0.f);
                    mu_sigma_stack[stack_index] = make_float4(0.f);
                    stack_index = stack_index-1;  //desapilo del stack.
                    //rtPrintf("", stack_index);

                    //if (launch_id == 128) rtPrintf("Saliendo del mesh - %d.\n", stack_index);
                    //if (launch_id == 128)
                    //rtPrintf("Stack index after sub: %d\n", stack_index);
                    //if (launch_id == 256) rtPrintf("Launch id: %d - Index: %d\n", launch_id, stack_index);
                } else {

                    if ((prd.mesh_index == attenuation_impedance_stack[stack_index-1].w) && (prd.mesh_index != EMPTY_STACK)) { //choco con el mesh que estaba afuera

                        attenuation_impedance_stack[stack_index-1] = attenuation_impedance_stack[stack_index]; //saco el material externo
                        mu_sigma_stack[stack_index-1] = mu_sigma_stack[stack_index];                           //y pongo el material que estaba
                        stack_index = stack_index-1;
                        //if (stack_index <= 0) rtPrintf("Estoy dentro de un vaso y estoy saliendo del material que lo rodea - %d.\n", stack_index);

                        //if (launch_id == 500) rtPrintf("Launch id: %d - Index: %d\n", launch_id, stack_index);

                        //if (launch_id == 128)
                        //rtPrintf("Stack index after sub: %d\n", stack_index);
                        //en el tope
                    } else { //choco con un nuevo mesh, estoy entrando.
                        if (attenuation_impedance_stack[stack_index].z){ //si estaba en un vaso
                            stack_index = stack_index+1;

                            //copio el vaso al tope
                            attenuation_impedance_stack[stack_index] = attenuation_impedance_stack[stack_index-1];
                            mu_sigma_stack[stack_index] = mu_sigma_stack[stack_index-1];

                            //agregar el nuevo abajo
                            attenuation_impedance_stack[stack_index-1].x = prd.media_attenuation;//prd.media_impedance,0.f,prd.mesh_index);
                            attenuation_impedance_stack[stack_index-1].y = prd.media_impedance;
                            attenuation_impedance_stack[stack_index-1].z = prd.is_vascular; //AGREGAR IS VASCULAR A LA PRD
                            attenuation_impedance_stack[stack_index-1].w = static_cast<float>(prd.mesh_index);
                            mu_sigma_stack[stack_index-1].x = prd.media_mu0;
                            mu_sigma_stack[stack_index-1].y = prd.media_mu1;
                            mu_sigma_stack[stack_index-1].z = prd.media_sigma;

                        } else {    //si no estaba en un vaso
                            stack_index = stack_index+1;  //apilo en el stack.
                            attenuation_impedance_stack[stack_index].x = prd.media_attenuation;//prd.media_impedance,0.f,prd.mesh_index);
                            attenuation_impedance_stack[stack_index].y = prd.media_impedance;
                            attenuation_impedance_stack[stack_index].w = static_cast<float>(prd.mesh_index);
                            mu_sigma_stack[stack_index].x = prd.media_mu0;// prd.media_mu1, prd.media_sigma, prd.mesh_index);
                            mu_sigma_stack[stack_index].y = prd.media_mu1;
                            mu_sigma_stack[stack_index].z = prd.media_sigma;
                        }
                    }
                }
            }

/*if (launch_id == 100){
                for (int i = 0; i <= 10; i++){
                    rtPrintf("%f - ", attenuation_impedance_stack[i].w);
                }
                rtPrintf("\n");}*/


            //if (launch_id == 128){
            //rtPrintf("stack_index: %d\n", stack_index);
            //rtPrintf("mesh_id: %d\n", attenuation_impedance_stack[stack_index].w);}
            //if (launch_id == 255) rtPrintf("stack handled\n");
            prd.depth++;
        }
        sample++;
        stack_index = 0;
        //if (launch_id == 256) rtPrintf("Launch id: %d - Index: %d\n", launch_id, stack_index);
        //if (launch_id == 512) rtPrintf("------ new frame ------");
    } while( sample < samples_te);

    /*//uint2 out_index;
    //uint3 dist_index;
    //uint3 imp_index;
    //rtPrintf("Distance: %f\n",prd.distance_traveled);
    //rtPrintf("Importance: %f\n", prd.importance);
    /*for (int i = 0; i < max_depth; i++){

        for (int sample = 0; sample < samples_te; sample++){
            dist_index = make_uint3(launch_id*max_depth+i, 0, sample);
            imp_index  = make_uint3(launch_id*max_depth+i, 1, sample);
            //if (results[imp_index] < intensity_epsilon)
            //    break;
            //if (launch_id == 255u)
            //    rtPrintf("Time elapsed: %f\n", results[dist_index]);
            float micros_traveled = results[dist_index];
            float row = micros_traveled / (axial_resolution / speed_of_sound);
            out_index  = make_uint2(launch_id, static_cast<unsigned int>(row)); //static_cast<unsigned int>(prd.distance_traveled));
            //rtPrintf("Micros_traveled: %f - Row: %f\n", micros_traveled, row);
            //rtPrintf("Index: %d, %d\n", out_index.x, out_index.y);
            //rtPrintf("Out_index: %d,%d - Imp_index: %d,%d - Importance: %f\n", out_index.x, out_index.y, imp_index.x, imp_index.y, results[imp_index]);
            if (out_index.y < max_rows){
                output_buffer[out_index] += make_float4(results[imp_index]);
            } else {
                output_buffer[make_uint2(launch_id, max_rows)] += make_float4(results[imp_index]);
            }

            results[imp_index] = 0;
        }

    }*/


}


RT_PROGRAM void exception(){
    const unsigned int code = rtGetExceptionCode();
    /*
    if ( code == RT_EXCEPTION_STACK_OVERFLOW )
        rtPrintf("La cague aca\n");

        //output_buffer[launch_id,static_cast<unsigned int>(prd.distance_traveled)] = error;
    else*/
    rtPrintExceptionDetails();
    //rtPrintExceptionDetails();


}


