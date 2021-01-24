#ifndef TRANSDUCER_H
#define TRANSDUCER_H

#include <units.h>
//#include <LinearMath/btVector3.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

#include <optixpp_namespace.h>
#include <optixu_math_stream_namespace.h>
#include <optixu_matrix_namespace.h>

#define M_PI 3.14159
#define X_AXIS optix::make_float3(1,0,0)
#define Y_AXIS optix::make_float3(0,1,0)
#define Z_AXIS optix::make_float3(0,0,1)

template<size_t transducer_elements>
class transducer                //light
{
public:
    struct transducer_element   //lights
    {
        optix::float3 position;
        optix::float3 direction;
    };

    transducer(const float frequency, const units::length::centimeter_t radius, units::length::millimeter_t transducer_element_separation,
               const optix::float3 & position, const std::array<units::angle::degree_t, 3> & angles) :
        frequency(frequency),
        radius(radius),
        position(position),
        angles(angles),
        transducer_element_separation(transducer_element_separation)
    {
        using namespace units::angle;
        using namespace units::literals;

        assert(transducer_element_separation * transducer_elements < M_PI * radius);

        radian_t x_angle { angles[0] }; //45º  ->  PI/4r
        radian_t y_angle { angles[1] }; //45º  ->  PI/4r
        radian_t z_angle { angles[2] }; //-90º -> -PI/2r

        auto amp = transducer_element_separation / radius;
        const radian_t amplitude { amp.to<float>() }; // angle covered by a single TE
        const radian_t angle_center_of_element { amplitude / 2.0f };

        radian_t angle = -(amplitude * transducer_elements / 2) + angle_center_of_element;

        for (size_t t = 0; t < transducer_elements; t++)
        {
            optix::float3 vector = optix::make_float3( std::sin( angle.to<float>()), std::cos( angle.to<float>()), 0 );
            rotate( vector, optix::Matrix4x4::rotate( z_angle.to<float>(), Z_AXIS ) );
            rotate( vector, optix::Matrix4x4::rotate( x_angle.to<float>(), X_AXIS ) );
            rotate( vector, optix::Matrix4x4::rotate( y_angle.to<float>(), Y_AXIS ) );
            elements[t] = transducer_element
            {
                position + radius.to<float>() * vector, // position
                vector  // direction
            };
            if (t == 0 || t == 511) std::cout << vector.x << "," << vector.y << "," << vector.z << std::endl;
            angle = angle + amplitude;
        }

    }

    void rotate( optix::float3& vector, optix::Matrix4x4 rotation_matrix )
    {
        float x, y, z;
        x = optix::dot( optix::make_float3( rotation_matrix.getRow(0) ), vector );
        y = optix::dot( optix::make_float3( rotation_matrix.getRow(1) ), vector );
        z = optix::dot( optix::make_float3( rotation_matrix.getRow(2) ), vector );
        vector.x = x; vector.y = y; vector.z = z;
    }

    transducer_element element(size_t i) const
    {
        return elements.at(i);
    }

    void print(bool direction) const
    {
        auto print_vec = [](const auto & v)
        {
            std::cout << v.x() << "," << v.z() << std::endl;
        };

        for (auto & element : elements)
        {
            print_vec(direction? element.direction : element.position);
        }
    }

    void update()
    {
        using namespace units::angle;
        using namespace units::literals;

        radian_t x_angle { angles[0] };
        radian_t y_angle { angles[1] };
        radian_t z_angle { angles[2] };

        auto amp = transducer_element_separation / radius;
        const radian_t amplitude { amp.template to<float>() }; // angle covered by a single TE
        const radian_t angle_center_of_element { amplitude / 2.0f };

        radian_t angle = -(amplitude * transducer_elements / 2) + angle_center_of_element;

        for (size_t t = 0; t < transducer_elements; t++)
        {
            optix::float3 vector = optix::make_float3( std::sin( angle.to<float>()), std::cos( angle.to<float>()), 0 );
            rotate( vector, optix::Matrix4x4::rotate( z_angle.to<float>(), Z_AXIS ) );
            rotate( vector, optix::Matrix4x4::rotate( x_angle.to<float>(), X_AXIS ) );
            rotate( vector, optix::Matrix4x4::rotate( y_angle.to<float>(), Y_AXIS ) );
            elements[t].position = position + radius.to<float>() * vector;

            elements[t].direction = vector;

            angle = angle + amplitude;
        }
    }

    void setPosition (const optix::float3 position)
    {
        this->position = position;
    }

    void setAngles(const std::array<units::angle::degree_t, 3> angles)
    {
        this->angles = angles;
    }

    optix::float3 getPosition ()
    {
        return position;
    }

    transducer_element &operator[](int index){
        return elements[index];
    }

    const float frequency;
    optix::float3 position, direction;
    const std::array<units::angle::degree_t, 3> angles;


private:
    const units::length::centimeter_t radius;
    const units::length::millimeter_t transducer_element_separation;


    std::array<transducer_element, transducer_elements> elements;
};

#endif // TRANSDUCER_H
