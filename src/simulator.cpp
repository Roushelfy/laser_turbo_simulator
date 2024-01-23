#include "simulator.h"
#include <sstream>
#include <fstream>
#include <filesystem>
#include <random>
Simulator::Simulator(const std::string &configFilePath, int argc, char *argv[]) : simulate_time(0), pwm(), laser(), pid(), record()
{
    // read params from config file
    std::ifstream configFile(configFilePath);
    std::string line;
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg.substr(0, 1) == "-")
        {
            auto delimiterPos = arg.find('=');
            auto key = arg.substr(1, delimiterPos - 1);
            auto value = arg.substr(delimiterPos + 1);
            config[key] = value;
        }
    }
    while (std::getline(configFile, line))
    {
        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, '='))
        {
            std::string value;
            std::getline(iss, value);
            if (key == "adc_frequency")
                adc_frequency = std::stod(value);
            else if (key == "galvo_delay")
                galvo_delay = std::stod(value);
            else if (key == "waist_radius")
                waist_radius = std::stod(value);
            else if (key == "waist_distance")
                waist_distance = std::stod(value);
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
            else if (key == "galvo_angular_speed")
                galvo_angular_speed = std::stod(value);
            else if (key == "fov")
                fov = std::stod(value);
            else if (key == "max_n_number")
                max_n_number = std::stoi(value);
            else if (key == "dac_resolution")
                dac_resolution = std::stoi(value);
            else if (key == "pwm_precision")
                pwm.precision = std::stod(value);
            else if (key == "long_time")
            {
                if (value == "true")
                    long_time = true;
                else
                    long_time = false;
            }
            else if (key == "debug")
            {
                if (value == "true")
                    debug = true;
                else
                    debug = false;
            }
            else if (key == "tracking_method")
            {
                if (value == "one_max_value")
                    method = tracking_method::one_max_value;
                else if (value == "n_max_value")
                    method = tracking_method::n_max_value;
                else if (value == "move_in_four")
                    method = tracking_method::move_in_four;
                else if (value == "fitting_gaussian")
                    method = tracking_method::fitting_gaussian;
                else
                    std::cout << "Tracking method not supported!" << std::endl;
            }
            else if (key == "pd_arrangement")
            {
                if (value == "square")
                    arrangement = pd_arrangement::square;
                else if (value == "circle")
                    arrangement = pd_arrangement::circle;
                else if (value == "cross")
                    arrangement = pd_arrangement::cross;
                else
                    std::cout << "PD arrangement not supported!" << std::endl;
            }
            else if (key == "object_radius")
            {
                object_radius = std::stod(value);
            }
            else if (key == "kp")
            {
                pid.kp = std::stod(value);
            }
            else if (key == "ki")
            {
                pid.ki = std::stod(value);
            }
            else if (key == "kd")
            {
                pid.kd = std::stod(value);
            }
            else if (key == "threshold")
            {
                pid.threshold = std::stod(value);
            }
            else if (key == "use_pid")
            {
                if (value == "true")
                    use_pid = true;
                else
                    use_pid = false;
            }
            else if (key == "use_log")
            {
                if (value == "true")
                    use_log = true;
                else
                    use_log = false;
            }
            else
                std::cout << "Key not supported!" << key << std::endl;
        }
    }
    for (auto &pair : config)
    {
        if (pair.first == "adc_frequency")
            adc_frequency = std::stod(pair.second);
        else if (pair.first == "galvo_delay")
            galvo_delay = std::stod(pair.second);
        else if (pair.first == "waist_radius")
            waist_radius = std::stod(pair.second);
        else if (pair.first == "waist_distance")
            waist_distance = std::stod(pair.second);
        else if (pair.first == "wavelength")
            wavelength = std::stod(pair.second);
        else if (pair.first == "central_intensity")
            central_intensity = std::stod(pair.second);
        else if (pair.first == "object_speed")
            object_speed = std::stod(pair.second);
        else if (pair.first == "object_moving_radius")
            object_moving_radius = std::stod(pair.second);
        else if (pair.first == "object_distance")
            object_distance = std::stod(pair.second);
        else if (pair.first == "pwm_precision")
            pwm.precision = std::stod(pair.second);
        else if (pair.first == "pd_number")
            pd_number = std::stod(pair.second);
        else if (pair.first == "pd_radius")
            pd_radius = std::stod(pair.second);
        else if (pair.first == "laser_position_x")
            laser_position(0) = std::stod(pair.second);
        else if (pair.first == "laser_position_y")
            laser_position(1) = std::stod(pair.second);
        else if (pair.first == "laser_position_z")
            laser_position(2) = std::stod(pair.second);
        else if (pair.first == "noise_stddev")
            noise_stddev = std::stod(pair.second);
        else if (pair.first == "background_intensity")
            background_intensity = std::stod(pair.second);
        else if (pair.first == "fov")
            fov = std::stod(pair.second);
        else if (pair.first == "max_n_number")
            max_n_number = std::stoi(pair.second);
        else if (pair.first == "dac_resolution")
            dac_resolution = std::stoi(pair.second);
        else if (pair.first == "galvo_angular_speed")
            galvo_angular_speed = std::stod(pair.second);
        else if (pair.first == "long_time")
        {
            if (pair.second == "true")
                long_time = true;
            else
                long_time = false;
        }
        else if (pair.first == "debug")
        {
            if (pair.second == "true")
                debug = true;
            else
                debug = false;
        }
        else if (pair.first == "kp")
        {
            pid.kp = std::stod(pair.second);
        }
        else if (pair.first == "ki")
        {
            pid.ki = std::stod(pair.second);
        }
        else if (pair.first == "kd")
        {
            pid.kd = std::stod(pair.second);
        }
        else if (pair.first == "threshold")
        {
            pid.threshold = std::stod(pair.second);
        }
        else if (pair.first == "use_pid")
        {
            if (pair.second == "true")
                use_pid = true;
            else
                use_pid = false;
        }
        else if (pair.first == "use_log")
        {
            if (pair.second == "true")
                use_log = true;
            else
                use_log = false;
        }
        else if (pair.first == "tracking_method")
        {
            if (pair.second == "one_max_value")
                method = tracking_method::one_max_value;
            else if (pair.second == "n_max_value")
                method = tracking_method::n_max_value;
            else if (pair.second == "move_in_four")
                method = tracking_method::move_in_four;
            else if (pair.second == "fitting_gaussian")
                method = tracking_method::fitting_gaussian;
            else
                std::cout << "Tracking method not supported!" << std::endl;
        }
        else if (pair.first == "pd_arrangement")
        {
            if (pair.second == "square")
                arrangement = pd_arrangement::square;
            else if (pair.second == "circle")
                arrangement = pd_arrangement::circle;
            else if (pair.second == "cross")
                arrangement = pd_arrangement::cross;
            else
                std::cout << "PD arrangement not supported!" << std::endl;
        }
        else if (pair.first == "object_radius")
        {
            object_radius = std::stod(pair.second);
        }
        else
            std::cout << "Key not supported!" << std::endl;
    }
    adc_time_step = 1.0 / adc_frequency;
    record = Record();
    object_angular_speed = object_speed / object_moving_radius;
    pwm.len = 1 / pwm.precision;
    // set laser to point to tag
    update_tag_position(0);
    laser.last_orientation = laser.target_orientation = (tag_position - laser_position).normalized();
    laser.last_time = laser.next_time = simulate_time;

    // set voltage based on orientation
    set_galvo_voltage();

    // generate pd data
    if (arrangement == pd_arrangement::square)
    {
        for (int i = 0; i < pd_number; ++i)
        {
            pd pd_temp;
            if (i < pd_number / 4)
            {
                pd_temp.position(0) = pd_radius;
                pd_temp.position(1) = -pd_radius + 2 * pd_radius * i / (pd_number / 4);
            }
            else if (i < pd_number / 2)
            {
                pd_temp.position(0) = pd_radius - 2 * pd_radius * (i - pd_number / 4) / (pd_number / 4);
                pd_temp.position(1) = pd_radius;
            }
            else if (i < 3 * pd_number / 4)
            {
                pd_temp.position(0) = -pd_radius;
                pd_temp.position(1) = pd_radius - 2 * pd_radius * (i - pd_number / 2) / (pd_number / 4);
            }
            else
            {
                pd_temp.position(0) = -pd_radius + 2 * pd_radius * (i - 3 * pd_number / 4) / (pd_number / 4);
                pd_temp.position(1) = -pd_radius;
            }
            pd_temp.value = 0;
            pd_temp.position(2) = 0;
            pd_array.push_back(pd_temp);
        }
        // print pd position
        if (debug == true)
            for (int i = 0; i < pd_number; ++i)
            {
                std::cout << "pd position:" << std::endl;
                std::cout << pd_array[i].position.transpose() << std::endl;
            }
    }
    else if (arrangement == pd_arrangement::cross)
    {
        double dr = pd_radius / (pd_number / 4);
        double r = dr;
        for (int i = 0; i < pd_number / 4; ++i)
        {
            pd pd_temp;
            pd_temp.position(0) = r;
            pd_temp.position(1) = 0;
            pd_temp.position(2) = 0;
            pd_temp.value = 0;
            pd_array.push_back(pd_temp);
            pd_temp.position(0) = -r;
            pd_array.push_back(pd_temp);
            pd_temp.position(0) = 0;
            pd_temp.position(1) = r;
            pd_array.push_back(pd_temp);
            pd_temp.position(1) = -r;
            pd_array.push_back(pd_temp);
            r += dr;
        }
    }
    else if (arrangement == pd_arrangement::circle)
    {
        for (int i = 0; i < pd_number; i++)
        {
            pd pd_temp;
            pd_temp.position(0) = pd_radius * std::cos(2 * M_PI * i / pd_number);
            pd_temp.position(1) = pd_radius * std::sin(2 * M_PI * i / pd_number);
            pd_temp.position(2) = 0;
            pd_temp.value = 0;
            pd_array.push_back(pd_temp);
        }
    }
    else
    {
        std::cout << "PD arrangement not supported!" << std::endl;
    }
    if (debug == true)
    {
        std::cout << "galvo_voltage:" << laser.galvo_voltage_x << " " << laser.galvo_voltage_y << std::endl;
        std::cout
            << "laser_orientation" << get_laser_orientation() << std::endl;
        std::cout << "tag_position" << tag_position.transpose() << std::endl;
        std::cout << "init dis" << distanceToLine(tag_position, laser_position, get_laser_orientation()) << std::endl;
        std::cout << "ADC frequency: " << adc_frequency << std::endl;
        std::cout << "ADC time step: " << adc_time_step << std::endl;
        std::cout << "Galvo delay: " << galvo_delay << std::endl;
        std::cout << "Waist radius: " << waist_radius << std::endl;
        std::cout << "Waist distance: " << waist_distance << std::endl;
        std::cout << "Wavelength: " << wavelength << std::endl;
        std::cout << "Central intensity: " << central_intensity << std::endl;
        std::cout << "Object speed: " << object_speed << std::endl;
        std::cout << "Object moving radius: " << object_moving_radius << std::endl;
        std::cout << "Object angular speed: " << object_angular_speed << std::endl;
        std::cout << "Object distance: " << object_distance << std::endl;
        std::cout << "PD number: " << pd_number << std::endl;
        std::cout << "PD radius: " << pd_radius << std::endl;
        std::cout << "Laser position: " << laser_position.transpose() << std::endl;
        std::cout << "Noise stddev: " << noise_stddev << std::endl;
        std::cout << "Background intensity: " << background_intensity << std::endl;
        std::cout << "FOV: " << fov << std::endl;
        std::cout << "Max n number: " << max_n_number << std::endl;
        std::cout << "Galvo angular speed: " << galvo_angular_speed << std::endl;
        std::cout << "DAC resolution: " << dac_resolution << std::endl;
        std::cout << "Long time: " << long_time << std::endl;
        std::cout << "Debug: " << debug << std::endl;
        std::cout << "use_log: " << use_log << std::endl;
        std::cout << "Use PID: " << use_pid << std::endl;

        std::cout << "PID kp: " << pid.kp << std::endl;
        std::cout << "PID ki: " << pid.ki << std::endl;
        std::cout << "PID kd: " << pid.kd << std::endl;
        std::cout << "PID threshold: " << pid.threshold << std::endl;
        std::cout << "Tracking method: ";
        switch (method)
        {
        case tracking_method::one_max_value:
            std::cout << "one_max_value" << std::endl;
            break;
        case tracking_method::n_max_value:
            std::cout << "n_max_value" << std::endl;
            break;
        case tracking_method::move_in_four:
            std::cout << "move_in_four" << std::endl;
            break;
        case tracking_method::fitting_gaussian:
            std::cout << "fitting_gaussian" << std::endl;
            break;
        default:
            break;
        }
        std::cout << "PD arrangement: ";
        switch (arrangement)
        {
        case pd_arrangement::square:
            std::cout << "square" << std::endl;
            break;
        case pd_arrangement::circle:
            std::cout << "circle" << std::endl;
            break;
        case pd_arrangement::cross:
            std::cout << "cross" << std::endl;
            break;
        default:
            break;
        }
        std::cout << "Object radius: " << object_radius << std::endl;
        std::cout << std::endl;
    }
}
Simulator::~Simulator()
{
}

void Simulator::run()
{
    if (debug == true)
        std::cout << "simulating..." << std::endl;
    while (true)
    {
        sample_pd();

        adjust_laser_orientation();

        // report every second
        if (std::fmod(simulate_time, 1) < adc_time_step * (pd_number) && simulate_time > 1)
        {
            if (debug == true)
            {
                std::cout << "Time: " << simulate_time << "s" << std::endl;
                std::cout << "Laser orientation: " << get_laser_orientation().transpose() << std::endl;
                std::cout << "Tag position: " << tag_position.transpose() << std::endl;
                std::cout << std::endl;
            }
            if (!long_time && simulate_time * object_angular_speed > 2 * M_PI)
            {
                if (use_log)
                {
                    // std::cout << "Tracking success!" << std::endl;
                    //  save record
                    // save all record data to files
                    std::string datapath = "../data/";
                    std::string distance = "dis_" + std::to_string((int)object_distance) + "/";
                    std::string speed = "speed_" + std::to_string((int)object_speed);
                    std::string filename_prefix = datapath + distance + speed;
                    std::filesystem::path path(filename_prefix);
                    if (!std::filesystem::exists(path))
                    {
                        // create
                        std::filesystem::create_directories(path);
                    }

                    std::string offsetfile = filename_prefix + "/offset.txt";
                    std::string snrfile = filename_prefix + "/snr.txt";
                    std::string recordfile = filename_prefix + "/record.txt";
                    std::ofstream recordFile(recordfile);
                    std::ofstream offsetFile(offsetfile);
                    std::ofstream snrFile(snrfile);
                    recordFile << "average_offset:" << record.average_distance << std::endl;
                    recordFile << "average_intensity:" << record.average_intensity << std::endl;
                    recordFile << "record_number:" << record.record_number << std::endl;
                    recordFile << "laser_offset:" << std::endl;
                    // 从最大值，最小值之间划分n个区间，统计次数
                    double max_dis = *record.laser_distance.rbegin();
                    double min_dis = *record.laser_distance.begin();
                    int plot_num = 100;
                    double interval = (max_dis - min_dis) / plot_num;
                    std::vector<int> count(plot_num, 0);
                    for (auto &dis : record.laser_distance)
                    {
                        int index = (dis - min_dis) / interval;
                        count[index]++;
                        offsetFile << dis << std::endl;
                    }
                    // sum all counts
                    int sum = 0;
                    for (auto &c : count)
                    {
                        sum += c;
                    }
                    for (int i = 0; i < plot_num; i++)
                    {
                        recordFile << min_dis + interval * i << ":" << count[i] / (double)sum << std::endl;
                    }
                    recordFile << "tag_intensity:" << std::endl;
                    double max_intensity = *record.tag_intensity.rbegin();
                    double min_intensity = *record.tag_intensity.begin();
                    double intensity_interval = (max_intensity - min_intensity) / plot_num;
                    std::vector<int> intensity_count(plot_num, 0);
                    for (auto &intensity : record.tag_intensity)
                    {
                        int index = (intensity - min_intensity) / intensity_interval;
                        intensity_count[index]++;
                    }
                    // sum all counts
                    sum = 0;
                    for (auto &c : intensity_count)
                    {
                        sum += c;
                    }
                    for (int i = 0; i < plot_num; i++)
                    {
                        recordFile << min_intensity + intensity_interval * i << ":" << intensity_count[i] / (double)sum << std::endl;
                    }
                    recordFile << "tag_SNR:" << std::endl;
                    double max_SNR = *record.tag_SNR.rbegin();
                    double min_SNR = *record.tag_SNR.begin();
                    double snr_interval = (max_SNR - min_SNR) / plot_num;
                    std::vector<int> snr_count(plot_num, 0);
                    for (auto &snr : record.tag_SNR)
                    {
                        int index = (snr - min_SNR) / snr_interval;
                        snr_count[index]++;
                        snrFile << snr << std::endl;
                    }
                    // sum all counts
                    sum = 0;
                    for (auto &c : snr_count)
                    {
                        sum += c;
                    }
                    for (int i = 0; i < plot_num; i++)
                    {
                        recordFile << min_SNR + snr_interval * i << ":" << snr_count[i] / (double)sum << std::endl;
                    }
                    std::cout << "object_speed:" << object_speed << ' ';
                    std::cout << "average_offset:" << record.average_distance << ' ';
                    std::cout << "average_intensity:" << record.average_intensity << std::endl;
                    recordFile.close();
                    offsetFile.close();
                    snrFile.close();
                }

                break;
            }
        }
        if (connection_lost())
        {
            std::cout << "Connection lost at: " << simulate_time << " seconds" << std::endl;
            std::cout << "object_speed:" << object_speed << std::endl;
            std::cout << "laser orientation" << get_laser_orientation().transpose() << std::endl;
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
    double w = waist_radius * std::sqrt(1 + std::pow((distance - waist_distance) / (M_PI * std::pow(waist_radius, 2) / wavelength), 2));
    double intensity = std::exp(-2 * std::pow(r, 2) / std::pow(w, 2)) * central_intensity * waist_radius * waist_radius / w / w;
    return intensity;
}

// based on galvo voltage(0-dac_resolution) ->(-0.5fov~+0.5fov)
// voltage applied based on z-axis
void Simulator::set_laser_orientation()
{
    double x = ((double)laser.galvo_voltage_x - dac_resolution / 2) / (dac_resolution / 2) * fov;
    double y = ((double)laser.galvo_voltage_y - dac_resolution / 2) / (dac_resolution / 2) * fov;
    laser.last_orientation = get_laser_orientation();
    laser.target_orientation(0) = std::sin(x) * std::cos(y);
    laser.target_orientation(1) = std::sin(y);
    laser.target_orientation(2) = std::cos(x) * std::cos(y);
    laser.last_time = simulate_time + galvo_delay;
    double time = std::max(abs(laser.target_orientation(0) - laser.last_orientation(0)), abs(laser.target_orientation(1) - laser.last_orientation(1))) / galvo_angular_speed;
    laser.next_time = simulate_time + galvo_delay + time;
    if (debug)
    {
        std::cout << "laser.last_orientation:" << laser.last_orientation.transpose() << std::endl;
        std::cout << "laser.target_orientation:" << laser.target_orientation.transpose() << std::endl;
        std::cout << "laser.last_time:" << laser.last_time << std::endl;
        std::cout << "laser.next_time:" << laser.next_time << std::endl;
    }
}

void Simulator::set_galvo_voltage()
{
    laser.galvo_avg_voltage_x = (laser.target_orientation(0) / fov) * dac_resolution / 2 + dac_resolution / 2;
    laser.galvo_avg_voltage_y = (laser.target_orientation(1) / fov) * dac_resolution / 2 + dac_resolution / 2;
    laser.galvo_voltage_x = round(laser.galvo_avg_voltage_x);
    laser.galvo_voltage_y = round(laser.galvo_avg_voltage_y);
    laser.last_galvo_voltage_x = laser.galvo_voltage_x;
    laser.last_galvo_voltage_y = laser.galvo_voltage_y;
}
void Simulator::adjust_laser_orientation()
{
    Eigen::Vector2d center;
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
        laser.galvo_voltage_x += pd_array[max_index].position(0) / small_edge;
        laser.galvo_voltage_y += pd_array[max_index].position(1) / small_edge;
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
        if (normalizer < 1e-4)
            normalizer = std::max(abs(x), abs(y));
        assert(normalizer != 0);
        laser.galvo_voltage_x += x / normalizer;
        laser.galvo_voltage_y += y / normalizer;
        break;
    }
    case tracking_method::move_in_four:
    {
        int max_index = 0;
        for (int i = 0; i < pd_number; i++)
        {
            if (pd_array[i].value > pd_array[max_index].value)
                max_index = i;
        }
        if (debug == true)
        {
            std::cout << "max_index:" << max_index << std::endl;
            std::cout << pd_array[max_index].position(0) << " " << pd_array[max_index].position(1) << std::endl;
        }
        if (abs(pd_array[max_index].position(0)) > abs(pd_array[max_index].position(1)))
        {
            laser.galvo_avg_voltage_x += pd_array[max_index].position(0) / abs(pd_array[max_index].position(0)) * pwm.precision;
        }
        else
        {
            laser.galvo_avg_voltage_y += pd_array[max_index].position(1) / abs(pd_array[max_index].position(1)) * pwm.precision;
        }
        if (debug == true)
        {
            std::cout << "galvo_avg_voltage_x:" << laser.galvo_avg_voltage_x << std::endl;
            std::cout << "galvo_avg_voltage_y:" << laser.galvo_avg_voltage_y << std::endl;
        }
        manage_pwm_queue();
        break;
    }
    case tracking_method::fitting_gaussian:
        center = find_center();
        if (debug == true)
        {
            std::cout << "center:" << center.transpose() << std::endl;
        }
        if (use_pid)
        {
            auto update_voltage = pid.update(center);
            laser.galvo_voltage_x += update_voltage.first;
            laser.galvo_voltage_y += update_voltage.second;
        }
        else
        {
            if (center(0) > 0)
            {
                if (laser.last_center(0) <= center(0))
                {
                    laser.galvo_voltage_x += 1;
                }
            }
            else if (center(0) < 0)
            {
                if (laser.last_center(0) >= center(0))
                {
                    laser.galvo_voltage_x -= 1;
                }
            }
            if (center(1) > 0)
            {
                if (laser.last_center(1) <= center(1))
                {
                    laser.galvo_voltage_y += 1;
                }
            }
            else if (center(1) < 0)
            {
                if (laser.last_center(1) >= center(1))
                {
                    laser.galvo_voltage_y -= 1;
                }
            }
            laser.last_center = center;
        }
        break;
    default:
        break;
    }
    // limit voltage
    if (laser.galvo_voltage_x > dac_resolution)
        laser.galvo_voltage_x = dac_resolution;
    if (laser.galvo_voltage_x < 0)
        laser.galvo_voltage_x = 0;
    if (laser.galvo_voltage_y > dac_resolution)
        laser.galvo_voltage_y = dac_resolution;
    if (laser.galvo_voltage_y < 0)
        laser.galvo_voltage_y = 0;
    set_laser_orientation();
}

void Simulator::sample_pd()
{
    for (int i = 0; i < pd_number / 2; i++)
    {
        update_tag_position(adc_time_step);
        Eigen::Vector3d new_laser_start = 2 * tag_position - laser_position;
        Eigen::Vector3d new_laser_orientation = get_laser_orientation() * -1;
        double intensity1 = background_intensity + gaussianBeamIntensity(pd_array[i].position, new_laser_start, new_laser_orientation);
        intensity1 += noise(noise_stddev);
        pd_array[i].value = std::max(intensity1 * 65535 > 65535 ? 65535 : intensity1 * 65535, (double)0);
        double intensity2 = background_intensity + gaussianBeamIntensity(pd_array[i + pd_number / 2].position, new_laser_start, new_laser_orientation);
        intensity2 += noise(noise_stddev);
        pd_array[i + pd_number / 2].value = std::max(intensity2 * 65535 > 65535 ? 65535 : intensity2 * 65535, (double)0);
    }
    // debug pd
    if (debug == true)
    {
        for (int i = 0; i < pd_number; i++)
        {
            std::cout << pd_array[i].value << " ";
        }
        std ::cout << std::endl;
    }
}

double Simulator::noise(double stddev)
{
    static std::mt19937 generator(std::random_device{}());
    std::normal_distribution<double> dist(0, stddev);
    return dist(generator);
}
Eigen::Vector3d Simulator::get_laser_orientation()
{
    Eigen::Vector3d res;
    if (simulate_time >= laser.next_time)
    {
        res = laser.target_orientation;
    }
    else if (simulate_time <= laser.last_time)
    {
        res = laser.last_orientation;
    }
    else
    {
        double time = simulate_time - laser.last_time;
        double ratio = time / (laser.next_time - laser.last_time);
        res = laser.last_orientation + (laser.target_orientation - laser.last_orientation) * ratio;
    }
    return res;
}
bool Simulator::connection_lost()
{
    double dis = distanceToLine(tag_position, laser_position, get_laser_orientation());
    double intensity = gaussianBeamIntensity(tag_position, laser_position, get_laser_orientation());
    // record
    if (use_log)
    {
        record.record_number++;
        record.laser_distance.insert(dis);
        record.sum_distance += dis;
        record.average_distance = record.sum_distance / record.record_number;
        record.tag_intensity.insert(intensity);
        record.sum_intensity += intensity;
        record.average_intensity = record.sum_intensity / record.record_number;
        double snr = 10 * log10(intensity / background_intensity);
        record.tag_SNR.insert(snr);
        record.sum_SNR += snr;
        record.average_SNR = record.sum_SNR / record.record_number;
    }
    if (dis > object_radius)
    {
        std::cout << "Connection lost! dis:" << dis << std::endl;
        std::cout << "object position:" << tag_position.transpose() << std::endl;
        std::cout << "laser pointing at" << get_laser_orientation().transpose() * object_distance << std::endl;
        return true;
    }
    else
        return false;
}

Eigen::Vector2d Simulator::find_center()
{
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // Define matrix A
    static MatrixXd A(pd_number, 3);
    static MatrixXd At;
    static MatrixXd AtA;
    static MatrixXd AtA_inv;
    // only calculate once
    if (AtA_inv.size() == 0)
    {
        for (int i = 0; i < pd_number; ++i)
        {
            A(i, 0) = pd_array[i].position(0);
            A(i, 1) = pd_array[i].position(1);
            A(i, 2) = 1.0;
        }

        // Compute transpose of A and A' * A
        At = A.transpose();
        AtA = At * A;
        AtA_inv = AtA.inverse();
    }
    // Prepare vector zfitVec
    VectorXd zfitVec(pd_number);
    for (int i = 0; i < pd_number; ++i)
    {
        zfitVec(i) = log(std::max(pd_array[i].value, 1));
    }

    // Calculate A' * zfitVec
    VectorXd Atzfit = At * zfitVec;

    if (debug)
    {
        std::cout << "A:" << std::endl;
        std::cout << A << std::endl;
        std::cout << "At:" << std::endl;
        std::cout << At << std::endl;
        std::cout << "AtA:" << std::endl;
        std::cout << AtA << std::endl;
        std::cout << "zfitVec:" << std::endl;
        std::cout << zfitVec << std::endl;
        std::cout << "Atzfit:" << std::endl;
        std::cout << Atzfit << std::endl;
    }

    // Solve the linear system (AtA) * aa = Atzfit
    VectorXd aa = AtA_inv * Atzfit;

    // Calculate x0 and y0
    double x0 = aa(0);
    double y0 = aa(1);

    // Print the center
    if (debug == true)
        std::cout << "Center: (" << x0 << ", " << y0 << ")" << std::endl;

    return Eigen::Vector2d(x0, y0);
}

void Simulator::manage_pwm_queue()
{
    std::pair<uint32_t, uint32_t> pwm_pair = {floor(laser.galvo_avg_voltage_x), floor(laser.galvo_avg_voltage_y)};
    if (debug == true)
        std::cout << "pwm_pair:" << pwm_pair.first << " " << pwm_pair.second << std::endl;
    if (pwm.average.first < laser.galvo_avg_voltage_x)
        pwm_pair.first++;
    if (pwm.average.second < laser.galvo_avg_voltage_y)
        pwm_pair.second++;
    pwm.queue.push(pwm_pair);
    pwm.sum.first += pwm_pair.first;
    pwm.sum.second += pwm_pair.second;
    if (pwm.queue.size() > pwm.len)
    {
        pwm.sum.first -= pwm.queue.front().first;
        pwm.sum.second -= pwm.queue.front().second;
        pwm.queue.pop();
    }
    pwm.average.first = (double)pwm.sum.first / pwm.queue.size();
    pwm.average.second = (double)pwm.sum.second / pwm.queue.size();
    laser.galvo_voltage_x = pwm_pair.first;
    laser.galvo_voltage_y = pwm_pair.second;
}
double Simulator::test_object_distance()
{
    double sum_dis = 0;
    int test_number = 100;
    object_angular_speed = 0;
    object_moving_radius = 0;
    laser.galvo_voltage_x = laser.galvo_voltage_y = dac_resolution / 2;
    double delta_angle = 2 * fov / dac_resolution;
    set_laser_orientation();
    update_tag_position(1);
    sample_pd();
    Eigen::Vector2d last_center = find_center();
    for (int i = 0; i < test_number; i++)
    {
        if (i % 4 == 0)
            laser.galvo_voltage_x += 1;
        else if (i % 4 == 1)
            laser.galvo_voltage_x -= 1;
        else if (i % 4 == 2)
            laser.galvo_voltage_y += 1;
        else
            laser.galvo_voltage_y -= 1;
        set_laser_orientation();
        update_tag_position(1);
        sample_pd();
        Eigen::Vector2d center = find_center();
        sum_dis += (center - last_center).norm() / tan(delta_angle);
        last_center = center;
    }
    double average_dis = sum_dis / test_number / 2;
    if (debug)
    {
        std::cout << "actual dis:" << object_distance << std::endl;
        std::cout << "tested average_dis:" << average_dis << std::endl;
        std::cout << "error " << (average_dis - object_distance) / object_distance * 100 << "%" << std::endl;
    }
    return average_dis;
}

void Simulator::test_dis_performance()
{
    int run_time = 1000;
    // run multiple times and measure stddev
    std::vector<double> dis;
    for (int i = 0; i < run_time; i++)
    {
        dis.push_back(test_object_distance());
    }
    double sum = 0;
    for (auto &d : dis)
    {
        sum += d;
    }
    double average = sum / run_time;
    double stddev = 0;
    for (auto &d : dis)
    {
        stddev += (d - average) * (d - average);
    }
    stddev = std::sqrt(stddev / run_time);
    std::cout << "average:" << average << std::endl;
    std::cout << "stddev:" << stddev << std::endl;
}
