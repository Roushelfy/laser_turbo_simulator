#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
uint32_t pd_number;
struct pd
{
    Eigen::Vector2d position;
    uint32_t value;
};
std::vector<pd> pd_array;
double x_0, y_0, sigma_x, sigma_y;
uint32_t pd_max_value;
double pd_radius;

Eigen::Vector2d find_center()
{
    using Eigen::HouseholderQR;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // Define matrix A
    MatrixXd A(pd_number, 3);
    for (int i = 0; i < pd_number; ++i)
    {
        A(i, 0) = pd_array[i].position(0);
        A(i, 1) = pd_array[i].position(1);
        A(i, 2) = 1;
    }
    // print matrix A
    // std::cout << "A = \n"
    //           << A << std::endl;

    // Compute transpose of A and A' * A
    MatrixXd At = A.transpose();
    MatrixXd AtA = At * A;

    // print AtA
    // std::cout << "AtA = \n"
    //           << AtA << std::endl;

    // Prepare vector zfitVec
    VectorXd zfitVec(pd_number);
    for (int i = 0; i < pd_number; ++i)
    {
        zfitVec(i) = log(pd_array[i].value);
    }

    // Calculate A' * zfitVec
    VectorXd Atzfit = At * zfitVec;

    // Solve the linear system (AtA) * aa = Atzfit
    // VectorXd aa = AtA.colPivHouseholderQr().solve(Atzfit);
    // // MatrixXd AtA_inv = AtA.completeOrthogonalDecomposition().pseudoInverse();
    MatrixXd AtA_inv = AtA.inverse();
    // //  print AtA_inv
    // //  std::cout << "AtA_inv = \n"
    // //            << AtA_inv << std::endl;

    // // VectorXd aa = AtA
    VectorXd aa = AtA_inv * Atzfit;

    // Calculate x0 and y0
    double a = aa(0);
    double b = aa(1);
    double c = aa(2);

    double logI = log(pd_max_value);
    // calculate delta
    double delta = 4 * c * c - 8 * c * logI + 4 * logI * logI - 4 * pd_radius * pd_radius * (a * a + b * b);

    //  Calculate sigma2
    double sigma2 = (2 * logI - 2 * c - sqrt(delta)) / (2 * (a * a + b * b));

    std::cout << "sigma2:" << sigma2 << std::endl;

    double x0 = a * sigma2;
    double y0 = b * sigma2;

    return Eigen::Vector2d(x0, y0);
}
double noise(double stddev)
{
    static std::mt19937 generator(std::random_device{}());
    std::normal_distribution<double> dist(0, stddev);
    return dist(generator);
}
void generate_pd_data()
{
    // genarate pd data based on x_0, y_0, sigma_x, sigma_y
    for (int i = 0; i < pd_number; ++i)
    {
        pd pd_temp;
        pd_temp.position(0) = pd_radius * std::cos(2 * M_PI * i / pd_number);
        pd_temp.position(1) = pd_radius * std::sin(2 * M_PI * i / pd_number);
        pd_temp.value = pd_max_value * std::exp(-(pd_temp.position(0) - x_0) * (pd_temp.position(0) - x_0) / (2 * sigma_x * sigma_x) - (pd_temp.position(1) - y_0) * (pd_temp.position(1) - y_0) / (2 * sigma_y * sigma_y));
        // add noise
        // pd_temp.value += rand() % 100 - 50;
        pd_temp.value += noise(0.01 * pd_max_value);
        pd_array.push_back(pd_temp);
    }

    // print pd value
    for (int i = 0; i < pd_number; ++i)
    {
        std::cout << pd_array[i].value << std::endl;
    }
}

void generate_another_pd_data()
{
    // genarate pd data around a square (n/4 on each side)
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
        pd_temp.value = pd_max_value * std::exp(-(pd_temp.position(0) - x_0) * (pd_temp.position(0) - x_0) / (2 * sigma_x * sigma_x) - (pd_temp.position(1) - y_0) * (pd_temp.position(1) - y_0) / (2 * sigma_y * sigma_y));
        // add noise
        pd_temp.value += noise(0.001 * pd_max_value);
        pd_array.push_back(pd_temp);
    }

    // print pd value and position
    for (int i = 0; i < pd_number; ++i)
    {
        std::cout << pd_array[i].value << std::endl;
        std::cout << pd_array[i].position.transpose() << std::endl;
    }
}
int main()
{
    // read pd data from file
    std::ifstream pd_file;
    pd_file.open("../pd.txt");
    if (!pd_file.is_open())
    {
        std::cout << "pd file open failed" << std::endl;
        return -1;
    }
    pd_file >> pd_number >> pd_radius;
    std::cout << pd_number << ' ' << pd_radius << std::endl;
    pd_file >> x_0 >> y_0 >> sigma_x >> sigma_y;
    std::cout << x_0 << ' ' << y_0 << ' ' << sigma_x << ' ' << sigma_y << std::endl;
    pd_file >> pd_max_value;
    std::cout << pd_max_value << std::endl;
    // generate_another_pd_data();
    generate_pd_data();
    pd_file.close();
    auto center = find_center();
    std::cout << "Center: (" << center.transpose() << ")" << std::endl;
}