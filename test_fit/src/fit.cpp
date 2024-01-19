#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
uint32_t pd_number;
struct pd
{
    Eigen::Vector2d position;
    uint16_t value;
};
std::vector<pd> pd_array;
double x_0, y_0, sigma_x, sigma_y;
uint16_t pd_max_value;
double pd_radius;

Eigen::Vector2d find_center()
{
    using Eigen::HouseholderQR;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // Define matrix A
    MatrixXd A(pd_number, 5);
    for (int i = 0; i < pd_number; ++i)
    {
        A(i, 0) = pd_array[i].position(0) * pd_array[i].position(0);
        A(i, 1) = pd_array[i].position(1) * pd_array[i].position(1);
        A(i, 2) = pd_array[i].position(0);
        A(i, 3) = pd_array[i].position(1);
        A(i, 4) = 1;
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
    MatrixXd AtA_inv = AtA.completeOrthogonalDecomposition().pseudoInverse();
    // print AtA_inv
    // std::cout << "AtA_inv = \n"
    //           << AtA_inv << std::endl;

    VectorXd aa = AtA_inv * Atzfit;

    // Calculate x0 and y0
    double x0 = -aa(2) / (2 * aa(0));
    double y0 = -aa(3) / (2 * aa(1));

    return Eigen::Vector2d(x0, y0);
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
        pd_temp.value += rand() % 100 - 50;
        pd_array.push_back(pd_temp);
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
    generate_pd_data();
    pd_file.close();
    std::cout << "Center: (" << find_center().transpose() << ")" << std::endl;
}