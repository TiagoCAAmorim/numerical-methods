#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

// Function to read a CSV file and return a matrix of doubles
std::vector<std::vector<double>> readCSV(const std::string& filePath) {
    std::vector<std::vector<double>> matrix;

    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return matrix;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::istringstream ss(line);
        std::string token;

        while (std::getline(ss, token, ',')) {
            try {
                row.push_back(std::stod(token));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid argument: " << e.what() << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Out of range: " << e.what() << std::endl;
            }
        }

        matrix.push_back(row);
    }

    file.close();

    return matrix;
}

std::vector<std::vector<double>> SolveGauss(std::vector<std::vector<double>> k, std::vector<std::vector<double>> f){

}

void tests_1D(){
    std::string n;
    fs::path k_filename;
    fs::path f_filename;
    fs::path x_true_filename;
    fs::path x_est_filename;

    fs::path currentPath = fs::current_path();
    std::cout << currentPath.string() << std::endl;

    std::vector<std::vector<double>> k_matrix;
    std::vector<std::vector<double>> f_vector;
    std::vector<std::vector<double>> x_vector;
    std::vector<std::vector<double>> x_true_vector;
    std::vector<std::vector<double>> x_est_vector;

    std::vector<int> a = {3, 5, 9, 10}; //, 20, 40, 50, 100, 200};
    for (int i : a) {
        n = std::to_string(i);
        std::cout << "Number of cells: " << n << std::endl;

        k_filename = currentPath / ("matrix/k_" + n + "_1.csv");
        f_filename = currentPath / ("matrix/f_" + n + "_1.csv");
        x_true_filename = currentPath / ("matrix/x_" + n + "_1.csv");
        x_est_filename = currentPath / ("matrix/x_prev_" + n + "_1.csv");

        k_matrix = readCSV(k_filename.string());
        f_vector = readCSV(f_filename.string());
        x_true_vector = readCSV(x_true_filename.string());
        x_est_vector = readCSV(x_est_filename.string());

        std::cout << "  K matrix = " << k_matrix.size() << " x " << k_matrix[0].size() << std::endl;
        std::cout << "  F vector = " << f_vector.size() << " x " << f_vector[0].size() << std::endl;
        std::cout << "  X vector = " << x_true_vector.size() << " x " << x_true_vector[0].size() << std::endl;
        std::cout << "  X est vector = " << x_est_vector.size() << " x " << x_est_vector[0].size() << std::endl;


    }
}


int main(){
    #ifdef _WIN32
        system("cls");
    #elif __unix__
        system("clear");
    #else
        std::cout << "Running on an unknown system" << std::endl;
    #endif

    tests_1D();
}