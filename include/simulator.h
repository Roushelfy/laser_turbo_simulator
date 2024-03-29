#ifndef SIMULATOR_H_
#define SIMULATOR_H_
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <set>
#include <map>
#include <queue>
#include <cmath>

struct pd
{
    Eigen::Vector3d position;
    int value;
};

struct laser_pointing
{
    Eigen::Vector3d last_orientation;
    Eigen::Vector3d target_orientation;
    uint32_t galvo_voltage_x;
    uint32_t galvo_voltage_y;
    uint32_t last_galvo_voltage_x;
    uint32_t last_galvo_voltage_y;
    Eigen::Vector2d last_center;
    double galvo_avg_voltage_x;
    double galvo_avg_voltage_y;
    double last_time;
    double next_time;
};
struct PID
{
    double kp;
    double ki;
    double kd;
    double threshold;
    Eigen::Vector2d last_error;
    Eigen::Vector2d sum_error;
    PID(double kp = 0, double ki = 0, double kd = 0) : kp(kp), ki(ki), kd(kd), last_error(), sum_error(), threshold(0) {}
    std::pair<int, int> update(const Eigen::Vector2d &error)
    {
        sum_error += error;
        Eigen::Vector2d delta_error = error - last_error;
        last_error = error;
        double delta_voltage_x = kp * error(0) + ki * sum_error(0) + kd * delta_error(0);
        double delta_voltage_y = kp * error(1) + ki * sum_error(1) + kd * delta_error(1);
        std::pair<int, int> delta_voltage = {0, 0};
        if (delta_voltage_x >= threshold)
            delta_voltage.first = 1;
        else if (delta_voltage_x <= -threshold)
            delta_voltage.first = -1;
        if (delta_voltage_y >= threshold)
            delta_voltage.second = 1;
        else if (delta_voltage_y <= -threshold)
            delta_voltage.second = -1;
        return delta_voltage;
    }
};
enum class pd_arrangement
{
    circle,
    square,
    cross,
};
enum class tracking_method
{
    one_max_value,
    n_max_value,
    move_in_four,
    fitting_gaussian,
};
struct pwm_queue
{
    std::queue<std::pair<uint32_t, uint32_t>> queue;
    std::pair<uint32_t, uint32_t> sum;
    std::pair<double, double> average;
    double precision;
    uint32_t len;
    pwm_queue() : sum({0, 0}), average({0, 0}) {}
};
struct Record
{
    std::set<double> laser_distance;
    double sum_distance;
    double average_distance;
    std::set<double> tag_intensity;
    std::set<double> tag_SNR;
    double sum_SNR;
    double average_SNR;
    double sum_intensity;
    double average_intensity;
    uint64_t record_number;
    Record() : average_distance(0), average_intensity(0), record_number(0), sum_distance(0), sum_intensity(0), sum_SNR(0), average_SNR(0) {}
};
class Simulator
{
private:
    /* data */
    std::map<std::string, std::string> config;
    double simulate_time;
    Eigen::Vector3d tag_position;
    Eigen::Vector3d laser_position;
    laser_pointing laser;
    tracking_method method;
    pd_arrangement arrangement;
    PID pid;
    pwm_queue pwm;
    double fov;
    double adc_frequency;  // 1/adc_frequency = adc time step
    double adc_time_step;  // sampling time between two adc
    double galvo_delay;    // time delay between adc and galvo
    double waist_radius;   // waist radius of laser in m
    double waist_distance; // waist distance of laser in m
    double wavelength;     // wavelength of laser in m
    double central_intensity;
    double object_speed;
    double object_moving_radius;
    double object_angular_speed;
    double object_distance;
    uint32_t pd_number;
    double pd_radius;
    double noise_stddev;
    double background_intensity;
    double object_radius;
    double galvo_angular_speed;
    bool debug;
    bool use_log;
    bool use_pid;
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
    void manage_pwm_queue();
    Eigen::Vector2d find_center();
    Eigen::Vector3d get_laser_orientation();
    double test_object_distance();
    void run();
    void test_dis_performance();
};
#endif