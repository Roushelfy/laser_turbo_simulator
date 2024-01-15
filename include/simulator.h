#ifndef SIMULATOR_H_
#define SIMULATOR_H_
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <set>
#include <map>
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
    move_in_four,
    fitting_galss,
};
struct Record
{
    std::set<double> laser_distance;
    double sum_distance;
    double average_distance;
    std::set<double> tag_intensity;
    double sum_intensity;
    double average_intensity;
    uint64_t record_number;
    Record() : average_distance(0), average_intensity(0), record_number(0), sum_distance(0), sum_intensity(0) {}
};
class Simulator
{
private:
    /* data */
    std::map<std::string, std::string> config;
    double simulate_time;
    Eigen::Vector3d laser_orientation;
    Eigen::Vector3d tag_position;
    Eigen::Vector3d laser_position;
    uint32_t galvo_voltage_x;
    uint32_t galvo_voltage_y;
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
    uint32_t dac_resolution;
    Record record;
    bool long_time;
    int max_n_number;
    std::vector<pd>
        pd_array;

public:
    Simulator(const std::string &configFilePath, int argc, char *argv[]);
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