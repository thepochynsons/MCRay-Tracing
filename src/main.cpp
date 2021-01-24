//#include "scene.h"
#include "volume.h"
#include "psf.h"
#include "rfimage.h"
#include "transducer.h"
#include "camera.h"
#include "random.cu"

#include <cmath>
#include <iostream>
#include <fstream>
#include <units.h>
#include <json.hpp>

//#include "inputmanager.h"
#include <memory>
#include <chrono>
#include <vector>
#include <random>
#include <algorithm>
#include <functional>
#include <nvrtc.h>
#include <optixpp_namespace.h>
#include <optixu_math_stream_namespace.h>

//
#define STRINGIFY(x) STRINGIFY2(x)
#define STRINGIFY2(x) #x
#define LINE_STR STRINGIFY(__LINE__)

// Error check/report helper for users of the C API
#define NVRTC_CHECK_ERROR( func )                                  \
  do {                                                             \
    nvrtcResult code = func;                                       \
    if( code != NVRTC_SUCCESS )                                    \
      throw optix::Exception( "ERROR: " __FILE__ "(" LINE_STR "): " +     \
          std::string( nvrtcGetErrorString( code ) ) );            \
  } while( 0 )

/*---Set units namespaces---*/
using namespace units::literals;
using namespace units::velocity;
using namespace units::length;
using namespace units::time;
using namespace units::angle;

/*-----Physics Settings-----*/
constexpr meters_per_second_t speed_of_sound = 1500_m / 1_s; // [μm/μs], [m/s]
constexpr float transducer_frequency = 4.5f; // [Mhz]
constexpr millimeter_t axial_resolution = millimeter_t(1.45f / transducer_frequency); // [mm], the division can be deduced from Burger13
constexpr size_t transducer_elements = 512;
constexpr size_t samples_te = 10; //santi
constexpr radian_t transducer_amplitude = 60_deg;
constexpr centimeter_t transducer_radius = 3_cm;
constexpr centimeter_t ultrasound_depth = 15_cm; // [1 5cm -> μm]
constexpr microsecond_t max_travel_time = microsecond_t(ultrasound_depth / speed_of_sound); // [μs]
constexpr unsigned int max_depth = 10u;  //Max number of ray bounces
unsigned int max_rows = (speed_of_sound * max_travel_time) / axial_resolution;
constexpr unsigned int resolution = 145; // [μm], from Burger13

/*-Static Types Definitions-*/
using psf_ = psf<7, 13, 7, resolution>;
using volume_ = volume<256, resolution>;
using rf_image_ = rf_image<transducer_elements, max_travel_time.to<unsigned int>(), static_cast<unsigned int>(axial_resolution.to<float>()*1000.0f/*mm->μm*/)>;
using transducer_ = transducer<transducer_elements>;



optix::Context context;
optix::Buffer normal_dist, uniform_dist;

void generate_random_numbers(std::array<unsigned int, transducer_elements> & seeds){
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<unsigned int> dist(0,UINT_MAX);
    for (unsigned int i=0; i < seeds.size(); i++){
        seeds[i] = dist(generator);
    }
}

void createContext(transducer_ transducer){
    context = optix::Context::create();
    int rtxon = 0;
    rtGlobalSetAttribute(RT_GLOBAL_ATTRIBUTE_ENABLE_RTX, sizeof (rtxon), &rtxon );
    context->setRayTypeCount(1);
    context->setEntryPointCount(1);

    context["samples_te"]->setUint(samples_te);
    context["max_rows"]->setUint(max_rows);
    context["transducer_frequency"]->setFloat(transducer_frequency);
    context["transducer_elements"]->setUint(transducer_elements);
    context["speed_of_sound"]->setFloat(speed_of_sound.to<float>());
    context["axial_resolution"]->setFloat(units::length::micrometer_t(axial_resolution).to<float>());
    context["max_travel_time"]->setFloat(max_travel_time.to<float>());

    //Lo que viene a continuacion es una bizarreada, tiene que poder hacerse mejor
    optix::float3 te_pos[transducer_elements];
    optix::float3 te_dir[transducer_elements];
    for (int i = 0; i < transducer_elements; i++){
        te_pos[i] = transducer[i].position;
        te_dir[i] = transducer[i].direction;
    }
    optix::Buffer transducer_positions_buffer = context->createBuffer(RT_BUFFER_INPUT);
    optix::Buffer transducer_directions_buffer = context->createBuffer(RT_BUFFER_INPUT);
    transducer_positions_buffer->setFormat( RT_FORMAT_USER );
    transducer_positions_buffer->setElementSize( sizeof( optix::float3 ) );
    transducer_positions_buffer->setSize( transducer_elements );
    transducer_directions_buffer->setFormat( RT_FORMAT_USER );
    transducer_directions_buffer->setElementSize( sizeof( optix::float3 ) );
    transducer_directions_buffer->setSize( transducer_elements );
    memcpy(transducer_positions_buffer->map(), &te_pos, sizeof( te_pos ));
    memcpy(transducer_directions_buffer->map(), &te_dir, sizeof( te_dir ));
    transducer_positions_buffer->unmap();
    transducer_directions_buffer->unmap();
    context["transducer_positions_buffer"]->set(transducer_positions_buffer);
    context["transducer_directions_buffer"]->set(transducer_directions_buffer);
    optix::Program ray_gen = context->createProgramFromPTXFile("../../ptx/ray_gen.ptx", "raygen" );
    context->setRayGenerationProgram( 0, ray_gen );

    optix::Program exception = context->createProgramFromPTXFile("../../ptx/ray_gen.ptx", "exception");
    context->setExceptionProgram( 0, exception );

    std::cout << "STACK SIZE: " <<  context->getStackSize() << std::endl;

    context->setExceptionProgram( 0, exception );
    context->setMissProgram(0, context->createProgramFromPTXFile("../../ptx/constantbg.ptx", "miss" ) );
    context["bg_color"]->setFloat(0.f);
    //context->setPrintEnabled( true );
    //context->setPrintBufferSize( 1048576 );

//-------------------------- Random values generation -----------------------------//
    //std::chrono::time_point<std::chrono::high_resolution_clock>startTime = std::chrono::high_resolution_clock::now();

    std::random_device rd1;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator1(rd1()); //Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<float> distribution1(0.0,1.0);
    std::array<float,transducer_elements*max_depth*samples_te> normal_distribution;
    std::generate(std::begin(normal_distribution), std::end(normal_distribution), std::bind(distribution1, generator1));
    //std::random_shuffle(normal_distribution.begin(), normal_distribution.end());

    std::random_device rd2;
    std::mt19937 generator2(rd2());
    std::uniform_real_distribution<float> distribution2(0.0,1.0);
    std::array<float,transducer_elements*max_depth*samples_te> uniform_distribution;
    std::generate_n(std::begin(uniform_distribution),transducer_elements*max_depth*samples_te,std::bind(distribution2,generator2));
    //std::generate(std::begin(uniform_distribution), std::end(uniform_distribution), std::bind(distribution2, generator2));
    //std::cout << transducer_elements*max_depth*samples_te << std::endl;
    /*std::cout<<"RANDOMS:"<<std::endl;
    for (int i=0; i< transducer_elements*max_depth*samples_te; i++){
        std::cout << uniform_distribution[i]<<std::endl;
    }*/


    normal_dist = context->createBuffer( RT_BUFFER_INPUT_OUTPUT );
    normal_dist->setFormat( RT_FORMAT_FLOAT );
    normal_dist->setSize(transducer_elements*samples_te*max_depth);
    uniform_dist = context->createBuffer( RT_BUFFER_INPUT_OUTPUT );
    uniform_dist->setFormat( RT_FORMAT_FLOAT );
    uniform_dist->setSize(transducer_elements*samples_te*max_depth);
    memcpy(normal_dist->map(), &normal_distribution, sizeof( normal_distribution ));
    memcpy(uniform_dist->map(), &uniform_distribution, sizeof( uniform_distribution ));
    normal_dist->unmap();
    uniform_dist->unmap();
    context["normal_dist"]->set(normal_dist);
    context["uniform_dist"]->set(uniform_dist);
    //std::chrono::time_point<std::chrono::high_resolution_clock>endTime = std::chrono::high_resolution_clock::now();
    //float timeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>
    //        (endTime-startTime).count();
    //std::cout << "Random generation Time: " <<  timeElapsed << "ms" << std::endl;
//---------------------------------------------------------------------------------//

    context["max_depth"]->setUint(max_depth);
    context["transudcer_elements"]->setUint(transducer_elements);

    optix::Buffer results = context->createBuffer( RT_BUFFER_OUTPUT , RT_FORMAT_INT , transducer_elements);   // 2 por distancia e intensidad.

    context["results"]->set(results);
}

int main(int argc, char** argv)
{
    // std::unique_ptr<inputManager> inputManager1;


    if (argc != 2)
    {
        std::cout << "Incorrect argument list." << std::endl;
        return 0;
    }

    static const volume_ texture_volume;

    const psf_ psf { transducer_frequency, 0.05f, .2f, 0.1f };

    rf_image_ rf_image { transducer_radius, transducer_amplitude, samples_te };

    nlohmann::json json;
    {
        std::ifstream infile { argv[1] };

        json << infile;
    }

    const auto & t_pos = json.at("transducerPosition");
    millimeter_t transducer_element_separation = transducer_amplitude.to<float>() * transducer_radius / transducer_elements;
    const auto & t_dir = json.at("transducerAngles");
    std::array<units::angle::degree_t, 3> transducer_angles = {degree_t((float)t_dir[0]), degree_t((float)t_dir[1]), degree_t((float)t_dir[2])};
    std::printf("t_pos: %f, %f, %f\nt_dir: %f, %f, %f\n", (float)t_pos[0],(float)t_pos[1],(float)t_pos[2],
                (float)t_dir[0],(float)t_dir[1],(float)t_dir[2]);

    transducer_ transducer(transducer_frequency, transducer_radius, transducer_element_separation,
                           optix::make_float3((float)t_pos[0], (float)t_pos[1], (float)t_pos[2]), transducer_angles);
    //btVector3(-17.0, 1.2, 6.45) // liver
    //btVector3(-19.5, 1.2, -0.45) // kidney

    createContext(transducer);


    /***
     * COMPILE RAY_GEN.CU
     * nvcc -I ./cuda,/home/santiago/Descargas/NVIDIA-OptiX-SDK-6.0.0-linux64/include/optixu,/home/santiago/Descargas/NVIDIA-OptiX-SDK-6.0.0-linux64/include,/home/santiago/Proyectos/MCRay-Tracing/include/units -L /home/santiago/Descargas/NVIDIA-OptiX-SDK-6.0.0-linux64/lib64/liboptix.so,/home/santiago/Descargas/NVIDIA-OptiX-SDK-6.0.0-linux64/lib64/liboptix_prime.so,/home/santiago/Descargas/NVIDIA-OptiX-SDK-6.0.0-linux64/lib64/liboptixu.so --ptx --output-file=./ptx/ray_gen.ptx ./src/ray_gen.cu
     *
    ***/
    try
    {
        loadScene(json, context);

        // Create InputManager
        //inputManager1 = std::make_unique<inputManager>();
        //inputManager1->setTransducer(& transducer);
        std::chrono::time_point<std::chrono::high_resolution_clock> RTstartTime, RTendTime,RMstartTime, RMendTime,PPstartTime, PPendTime;

//
        /*----------Vars Init----------*/
        optix::Buffer output_buffer, segments_buffer;
        optix::float3 point;
        ray_physics::segment* buffer_base;
        float RTtimeElapsed, RMtimeElapsed,PPtimeElapsed,intensity, scattering;
        unsigned int ray_i, depth, i, steps, step, sample;
        units::time::microsecond_t time_elapsed;
        ray_physics::segment segment;

        //int iteracion = 0;

        while(true)
        {


            output_buffer = sutil::createOutputBuffer(context, RT_FORMAT_FLOAT4, transducer_elements, max_rows, false);
            context["output_buffer"]->set(output_buffer);
            segments_buffer = context->createBuffer(RT_BUFFER_OUTPUT, RT_FORMAT_USER, transducer_elements * max_depth * samples_te); //
            segments_buffer->setElementSize( sizeof(ray_physics::segment) );
            context["segments_buffer"]->set(segments_buffer);


            // Calculate time elapsed from last loop



            // Update Input
            //inputManager1->update(timeElapsed);

            rf_image.clear();
            RTstartTime = std::chrono::high_resolution_clock::now();
            context->launch( 0, transducer_elements  );
            RTendTime = std::chrono::high_resolution_clock::now();


            RMstartTime = std::chrono::high_resolution_clock::now();
            buffer_base= static_cast<ray_physics::segment*>( segments_buffer->map());
            for (sample = 0; sample < samples_te; sample++){
                for (ray_i = 0; ray_i < transducer_elements; ray_i++)
                {
                    for (depth = 0; depth < max_depth; depth++)
                    {

                        //i = (ray_i * max_depth + depth) * samples_te + sample;
                        i = ray_i + transducer_elements * (depth + max_depth * sample);
                        //std::printf("ray_i: %d - depth: %d - sample: %d - max_depth: %d - samples_te: %d --> Index: %d\n",ray_i,depth,sample,max_depth,samples_te,i);
                        segment = buffer_base[i/*(ray_i * max_depth + depth)+(sample*max_depth*transducer_elements)*/];
                        //if (!((segment.from.x == segment.to.x) && (segment.from.y == segment.to.y) && (segment.from.z == segment.to.z)))

                        const auto starting_micros = rf_image.micros_traveled(units::length::micrometer_t(segment.distance_traveled*1000.f)); //mm -> μm
                        const auto d = distance(segment.from, segment.to);//.to<float>(); // [mm]

                        steps = d / (axial_resolution );

                        const auto delta_step = /*units::length::micrometer_t(axial_resolution).to<float>()*/axial_resolution.to<float>() * segment.direction;
                        const auto time_step = rf_image.micros_traveled(axial_resolution); // [μs]
                        point = segment.from;
                        time_elapsed = starting_micros;
                        intensity = segment.initial_intensity;
                        //std::cout << "depth " << depth <<  " - segment initial intensity " <<segment.initial_intensity << std::endl;
                        for (step = 0; step < steps && time_elapsed < max_travel_time; step++)
                        {
                            scattering = texture_volume.get_scattering(segment.mu1, segment.mu0, segment.sigma, point.x, point.y, point.z);

                            rf_image.add_echo(ray_i, intensity * scattering, time_elapsed);

                            // Step forward through the segment, decreasing intensity using Beer-Lambert's law
                            point = point + delta_step;
                            time_elapsed = time_elapsed + time_step;

                            constexpr auto k = 1.0f;
                            intensity *= std::exp(-segment.attenuation * axial_resolution.to<float>() * 0.01f * transducer_frequency ); // *k

                        }
                        // Add reflection term, i.e. intensity directly reflected back to the transducer. See Burger13, Eq. 10.

                        rf_image.add_echo(ray_i, segment.reflected_intensity, starting_micros + time_step * (steps-1));
                        //if (segment.reflected_intensity) std::cout << segment.reflected_intensity << std::endl;
                        //std::cout << "ray_i: " << ray_i << " - depth: " << depth << " - sample: " << sample << " - BTTI: " << segment.reflected_intensity << std::endl;
                            //if (ray_i == 190) std::printf("Final echo\npoint: %f,%f,%f - time_elapsed: %f\n",segment.to.x, segment.to.y, segment.to.z,time_elapsed);
                            //std::cout << "Segment from " << segment.from << " to " << segment.to << std::endl;


                        //std::cout << "For sample" << std::endl;
                    }//}
                    //std::cout << "For depth" << std::endl;
                }
            }
            //std::cout << "For ray_i" << std::endl;
            segments_buffer->unmap();
            //std::cout << "echoes added: " << rf_image.getIterations() << std::endl;
            RMendTime = std::chrono::high_resolution_clock::now();


            PPstartTime = std::chrono::high_resolution_clock::now();

            rf_image.convolve(psf);
            rf_image.envelope();
            rf_image.postprocess();

            PPendTime = std::chrono::high_resolution_clock::now();

            RMtimeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                    (RMendTime-RMstartTime).count();
            RTtimeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                    (RTendTime-RTstartTime).count();
            PPtimeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                    (PPendTime-PPstartTime).count();

            printf("Times: %fms, %fms, %fms, %fms\n",RTtimeElapsed,RMtimeElapsed,PPtimeElapsed,RTtimeElapsed+RMtimeElapsed+PPtimeElapsed);
            int* results_base= static_cast<int*>( context["results"]->getBuffer()->map());
            int sum = 0;
            for (int i = 0; i<transducer_elements; i++){
                sum = sum + results_base[i];
            }
            context["results"]->getBuffer()->unmap();
            printf("Rayos: %d\n", sum);
            rf_image.show();
            context["output_buffer"]->getBuffer()->destroy();
            context["segments_buffer"]->getBuffer()->destroy();
            //std::cout << "Context launched" << std::endl;



//            auto rays = scene.cast_rays<samples_te,transducer_elements>(transducer);
//            for (unsigned int ray_i = 0; ray_i < rays.size(); ray_i++)
//            {
//                const auto & ray = rays[ray_i];
//                for (unsigned int sample_i = 0; sample_i < samples_te; sample_i++)
//                {
//                    const auto & sample = ray[sample_i];
//                    for (auto & segment : sample)
//                    {
//                        // Add reflection term, i.e. intensity directly reflected back to the transducer. See Burger13, Eq. 10.
//                        rf_image.add_echo(ray_i, (segment.reflected_intensity)/samples_te, starting_micros + time_step * (steps-1));
//                    }
//                }
//            }


            //std::cout << "ms - Ray marching+postprocessing time: " <<  timeElapsed << "ms - FPS: " << 1000.f/timeElapsed << std::endl;
            //std::cout << "Frame Time: " <<  timeElapsed << "ms - FPS: " << 1000.f/timeElapsed << std::endl;
        }
    }
    catch (const std::exception & ex)
    {
        std::cout << "The program found an error and will terminate.\n"
                  << "Reason:\n"
                  << ex.what() << std::endl;
    }




}
