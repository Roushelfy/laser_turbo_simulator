#ifndef SIMULATOR_H_
#define SIMULATOR_H_
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <cmath>

struct pd
{
    Eigen::Vector3d position;
    uint16_t value;
};

enum class tracking_method
{
    one_max_value,
    n_max_value,
    fitting_galss,
};
class Simulator
{
private:
    /* data */
    double simulate_time;
    Eigen::Vector3d laser_orientation;
    Eigen::Vector3d tag_position;
    Eigen::Vector3d laser_position;
    uint16_t galvo_voltage_x;
    uint16_t galvo_voltage_y;
    tracking_method method;
    double time_step;
    double fov;
    double adc_frequency; // 1/adc_frequency = adc time step
    double adc_time_step; // sampling time between two adc
    double galvo_delay;   // time delay between adc and galvo
    double waist_radius;  // waist radius of laser in m
    double wavelength;    // wavelength of laser in m
    double central_intensity;
    double object_speed;
    double object_moving_radius;
    double object_angular_speed;
    double object_distance;
    double pd_number;
    double pd_radius;
    double noise_stddev;
    double background_intensity;
    double object_radius;
    int max_n_number;
    std::vector<pd>
        pd_array;

public:
    Simulator(const std::string &configFilePath);
    ~Simulator();
    void set_laser_orientation();
    void set_galvo_voltage();
    void update_tag_position(double time_pass);
    void sample_pd();
    void adjust_laser_orientation();
    double distanceToLine(const Eigen::Vector3d &point,
                          const Eigen::Vector3d &pstart,
                          const Eigen::Vector3d &orientation);
    double gaussianBeamIntensity(const Eigen::Vector3d &dest,
                                 const Eigen::Vector3d &start,
                                 const Eigen::Vector3d &orientation);
    double noise(double stddev);
    bool connection_lost();
    void run();
};
#endif