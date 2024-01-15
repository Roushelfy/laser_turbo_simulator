#include "simulator.h"
#include <sstream>
#include <random>
Simulator::Simulator(const std::string &configFilePath) : simulate_time(0)
{
    // read params from config file
    std::ifstream configFile(configFilePath);
    std::string line;
    while (std::getline(configFile, line))
    {
        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, '='))
        {
            std::string value;
            std::getline(iss, value);
            if (key == "time_step")
                time_step = std::stod(value);
            else if (key == "adc_frequency")
                adc_frequency = std::stod(value);
            else if (key == "galvo_delay")
                galvo_delay = std::stod(value);
            else if (key == "waist_radius")
                waist_radius = std::stod(value);
            else if (key == "wavelength")
                wavelength = std::stod(value);
            else if (key == "central_intensity")
                central_intensity = std::stod(value);
            else if (key == "object_speed")
                object_speed = std::stod(value);
            else if (key == "object_moving_radius")
                object_moving_radius = std::stod(value);
            else if (key == "object_distance")
                object_distance = std::stod(value);
            else if (key == "pd_number")
                pd_number = std::stod(value);
            else if (key == "pd_radius")
                pd_radius = std::stod(value);
            else if (key == "laser_position_x")
                laser_position(0) = std::stod(value);
            else if (key == "laser_position_y")
                laser_position(1) = std::stod(value);
            else if (key == "laser_position_z")
                laser_position(2) = std::stod(value);
            else if (key == "noise_stddev")
                noise_stddev = std::stod(value);
            else if (key == "background_intensity")
                background_intensity = std::stod(value);
            else if (key == "fov")
                fov = std::stod(value);
            else if (key == "max_n_number")
                max_n_number = std::stoi(value);
            else if (key == "long_time")
            {
                if (value == "true")
                    long_time = true;
                else
                    long_time = false;
            }
            else if (key == "tracking_method")
            {
                if (value == "one_max_value")
                    method = tracking_method::one_max_value;
                else if (value == "three_max_value")
                    method = tracking_method::n_max_value;
                else if (value == "fitting_glass")
                    method = tracking_method::fitting_galss;
                else
                    std::cout << "Tracking method not supported!" << std::endl;
            }
            else if (key == "object_radius")
            {
                object_radius = std::stod(value);
            }
            else
                std::cout << "Key not supported!" << std::endl;
        }
    }
    adc_time_step = 1.0 / adc_frequency;
    record = Record();
    object_angular_speed = object_speed / object_moving_radius;
    // set laser to point to tag
    update_tag_position(0);
    laser_orientation = (tag_position - laser_position).normalized();
    // set voltage based on orientation
    set_galvo_voltage();
    std::cout
        << "laser_orientation" << laser_orientation.transpose() << std::endl;
    std::cout << "tag_position" << tag_position.transpose() << std::endl;
    std::cout << "init dis" << distanceToLine(tag_position, laser_position, laser_orientation) << std::endl;
    for (int i = 0; i < pd_number; i++)
    {
        pd pd_temp;
        pd_temp.position(0) = pd_radius * std::cos(2 * M_PI * i / pd_number);
        pd_temp.position(1) = pd_radius * std::sin(2 * M_PI * i / pd_number);
        pd_temp.position(2) = 0;
        pd_temp.value = 0;
        pd_array.push_back(pd_temp);
    }
    std::cout << "Time step: " << time_step << std::endl;
    std::cout << "ADC frequency: " << adc_frequency << std::endl;
    std::cout << "ADC time step: " << adc_time_step << std::endl;
    std::cout << "Galvo delay: " << galvo_delay << std::endl;
    std::cout << "Waist radius: " << waist_radius << std::endl;
    std::cout << "Wavelength: " << wavelength << std::endl;
    std::cout << "Central intensity: " << central_intensity << std::endl;
    std::cout << "Object speed: " << object_speed << std::endl;
    std::cout << "Object moving radius: " << object_moving_radius << std::endl;
    std::cout << "Object angular speed: " << object_angular_speed << std::endl;
    std::cout << "Object distance: " << object_distance << std::endl;
    std::cout << "PD number: " << pd_number << std::endl;
    std::cout << "PD radius: " << pd_radius << std::endl;
    std::cout << "Laser position: " << laser_position.transpose() << std::endl;
}
Simulator::~Simulator()
{
}

void Simulator::run()
{
    std::cout << "simulating..." << std::endl;
    while (true)
    {
        sample_pd();

        adjust_laser_orientation();

        // report every second
        if (std::fmod(simulate_time, 1) < 1e-5 && simulate_time > 1)
        {
            std::cout << "Time: " << simulate_time << "s" << std::endl;
            std::cout << "Laser orientation: " << laser_orientation.transpose() << std::endl;
            std::cout << "Tag position: " << tag_position.transpose() << std::endl;
            std::cout << std::endl;
            if (!long_time)
            {
                std::cout << "Tracking success!" << std::endl;
                // save record
                std::ofstream recordFile("../data/record.txt");
                recordFile << "average_distance:" << record.average_distance << std::endl;
                recordFile << "average_intensity:" << record.average_intensity << std::endl;
                recordFile << "record_number:" << record.record_number << std::endl;
                recordFile << "laser_distance:" << std::endl;
                // 从最大值，最小值之间划分20个区间，统计次数
                double max_dis = *record.laser_distance.rbegin();
                double min_dis = *record.laser_distance.begin();
                double interval = (max_dis - min_dis) / 20;
                std::vector<int> count(20, 0);
                for (auto &dis : record.laser_distance)
                {
                    int index = (dis - min_dis) / interval;
                    count[index]++;
                }
                for (int i = 0; i < 20; i++)
                {
                    recordFile << min_dis + interval * i << ":" << count[i] << std::endl;
                }
                recordFile << "tag_intensity:" << std::endl;
                double max_intensity = *record.tag_intensity.rbegin();
                double min_intensity = *record.tag_intensity.begin();
                double intensity_interval = (max_intensity - min_intensity) / 20;
                std::vector<int> intensity_count(20, 0);
                for (auto &intensity : record.tag_intensity)
                {
                    int index = (intensity - min_intensity) / intensity_interval;
                    intensity_count[index]++;
                }
                for (int i = 0; i < 20; i++)
                {
                    recordFile << min_intensity + intensity_interval * i << ":" << intensity_count[i] << std::endl;
                }
                recordFile.close();
                // 绘制分布曲线
                break;
            }
        }
        if (connection_lost())
        {
            std::cout << "Connection lost at: " << simulate_time << " seconds" << std::endl;
            std::cout << "laser orientation" << laser_orientation.transpose() << std::endl;
            std::cout << "tag position" << tag_position.transpose() << std::endl;
            break;
        }
    }
}

void Simulator::update_tag_position(double time_pass)
{
    simulate_time += time_pass;
    tag_position(0) = object_moving_radius * std::cos(object_angular_speed * simulate_time);
    tag_position(1) = object_moving_radius * std::sin(object_angular_speed * simulate_time);
    tag_position(2) = object_distance;
    // if (time_pass != 0)
    // std::cout << distanceToLine(tag_position, laser_position, laser_orientation) << std::endl;
}

double Simulator::distanceToLine(const Eigen::Vector3d &point,
                                 const Eigen::Vector3d &pstart,
                                 const Eigen::Vector3d &orientation)
{
    // make sure orientation is not zero
    assert(orientation.norm() != 0);
    Eigen::Vector3d diff = point - pstart;
    Eigen::Vector3d crossProd = diff.cross(orientation);
    return crossProd.norm() / orientation.norm();
}

double Simulator::gaussianBeamIntensity(const Eigen::Vector3d &dest,
                                        const Eigen::Vector3d &start,
                                        const Eigen::Vector3d &orientation)
{
    double distance = (dest - start).norm();
    double r = distanceToLine(dest, start, orientation);
    double w = waist_radius * std::sqrt(1 + std::pow(distance / (M_PI * std::pow(waist_radius, 2) / wavelength), 2));
    double intensity = std::exp(-2 * std::pow(r, 2) / std::pow(w, 2)) * central_intensity;
    return intensity;
}

// based on galvo voltage(0-4095) ->(-fov~+fov)
// voltage applied based on z-axis
void Simulator::set_laser_orientation()
{
    double x = (galvo_voltage_x - 2048) / 2048.0 * fov;
    double y = (galvo_voltage_y - 2048) / 2048.0 * fov;
    laser_orientation(0) = std::sin(x) * std::cos(y);
    laser_orientation(1) = std::sin(y);
    laser_orientation(2) = std::cos(x) * std::cos(y);
}

void Simulator::set_galvo_voltage()
{
    galvo_voltage_x = (laser_orientation(0) / fov) * 2048 + 2048;
    galvo_voltage_y = laser_orientation(1) / fov * 2048 + 2048;
}
void Simulator::adjust_laser_orientation()
{
    update_tag_position(galvo_delay);
    switch (method)
    {
    case tracking_method::one_max_value:
    {
        int max_index = 0;
        for (int i = 0; i < pd_number; i++)
        {
            if (pd_array[i].value > pd_array[max_index].value)
                max_index = i;
        }
        // normalize;
        double small_edge = std::min(abs(pd_array[max_index].position(0)), abs(pd_array[max_index].position(1)));
        if (small_edge < 1e-4)
            small_edge = std::max(abs(pd_array[max_index].position(0)), abs(pd_array[max_index].position(1)));
        assert(small_edge != 0);
        galvo_voltage_x += pd_array[max_index].position(0) / small_edge;
        galvo_voltage_y += pd_array[max_index].position(1) / small_edge;
        break;
    }
    case tracking_method::n_max_value:
    {
        std::vector<pd> temporate_pd;
        for (int i = 0; i < pd_number; i++)
        {
            temporate_pd.push_back(pd_array[i]);
        }
        std::sort(temporate_pd.begin(), temporate_pd.end(), [](const pd &a, const pd &b)
                  { return a.value > b.value; });
        double x = 0;
        double y = 0;
        for (int i = 0; i < max_n_number; i++)
        {
            x += temporate_pd[i].position(0) * temporate_pd[i].value;
            y += temporate_pd[i].position(1) * temporate_pd[i].value;
        }
        double normalizer = std::min(abs(x), abs(y));
        if (normalizer == 0)
            normalizer = std::max(abs(x), abs(y));
        assert(normalizer != 0);
        galvo_voltage_x += x / normalizer;
        galvo_voltage_y += y / normalizer;
        break;
    }
    case tracking_method::fitting_galss:
    default:
        break;
    }
    if (galvo_voltage_x > 4095)
        galvo_voltage_x = 4095;
    if (galvo_voltage_x < 0)
        galvo_voltage_x = 0;
    if (galvo_voltage_y > 4095)
        galvo_voltage_y = 4095;
    if (galvo_voltage_y < 0)
        galvo_voltage_y = 0;
    set_laser_orientation();
}

void Simulator::sample_pd()
{
    for (int i = 0; i < pd_number / 2; i++)
    {
        update_tag_position(adc_time_step);
        Eigen::Vector3d new_laser_start = (2 * tag_position - laser_position) * 0.5;
        Eigen::Vector3d new_laser_orientation = laser_orientation * -1;
        double intensity1 = background_intensity + gaussianBeamIntensity(pd_array[i].position, new_laser_start, new_laser_orientation);
        intensity1 += noise(noise_stddev);
        pd_array[i].value = intensity1 * 65535 > 65535 ? 65535 : intensity1 * 65535;
        double intensity2 = background_intensity + gaussianBeamIntensity(pd_array[i + pd_number / 2].position, new_laser_start, new_laser_orientation);
        intensity2 += noise(noise_stddev);
        pd_array[i + pd_number / 2].value = intensity2 * 65535 > 65535 ? 65535 : intensity2 * 65535;
    }
}

double Simulator::noise(double stddev)
{
    static std::mt19937 generator(std::random_device{}());
    std::normal_distribution<double> dist(0, stddev);
    return dist(generator);
}

bool Simulator::connection_lost()
{
    double dis = distanceToLine(tag_position, laser_position, laser_orientation);
    double intensity = gaussianBeamIntensity(tag_position, laser_position, laser_orientation);
    // record
    record.record_number++;
    record.laser_distance.insert(dis);
    record.sum_distance += dis;
    record.average_distance = record.sum_distance / record.record_number;
    record.tag_intensity.insert(intensity);
    record.sum_intensity += intensity;
    record.average_intensity = record.sum_intensity / record.record_number;
    if (dis > object_radius)
    {
        std::cout << "Connection lost! dis:" << dis << std::endl;
        std::cout << "object position:" << tag_position.transpose() << std::endl;
        return true;
    }
    else
        return false;
}