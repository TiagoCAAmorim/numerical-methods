#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>
#include <chrono>
#include <cmath>

// using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

const double eps = 1e-17;


void printMatrix(std::vector<std::vector<double>> matrix, std::ostream& output){
    for (size_t r=0; r<matrix.size(); r++){
        std::string sep = "| ";
        for (size_t c=0; c<matrix[r].size(); c++){
            output << sep << matrix[r][c];
            sep = ", ";
        }
        output << " |" << std::endl;
    }
}

void printVector(std::vector<double> vector, std::ostream& output){
    std::string sep = "{";
    for (size_t i=0; i<vector.size(); i++){
        output << sep << vector[i];
        sep = ", ";
    }
    output << "}" << std::endl;
}

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

std::vector<double> flattenMatrix(const std::vector<std::vector<double>>& matrix) {
    std::vector<double> flattenedVector;

    for (const auto& row : matrix) {
        flattenedVector.insert(flattenedVector.end(), row.begin(), row.end());
    }
    return flattenedVector;
}

std::vector<std::vector<double>> combineMatrixAndVector(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& vector) {

    if (matrix.size() != vector.size() || matrix.empty() || matrix[0].empty()) {
        throw std::invalid_argument("Error: Invalid matrix or vector dimensions");
    }

    size_t numRows = matrix.size();
    size_t numCols = matrix[0].size() + 1;
    std::vector<std::vector<double>> result(numRows, std::vector<double>(numCols));

    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            result[i][j] = matrix[i][j];
        }
        result[i][numCols - 1] = vector[i];
    }

    return result;
}

std::vector<double> matMult(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& vector) {

    size_t numRows = matrix.size();
    size_t numCols = matrix[0].size();
    if (numCols != vector.size()) {
        throw std::invalid_argument("Error: Matrix and vector dimensions are incompatible for multiplication");
    }

    std::vector<double> result(numRows, 0.0);
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

std::vector<double> vecDif(
    const std::vector<double>& vector1,
    const std::vector<double>& vector2) {

    if (vector1.size() != vector2.size()) {
        throw std::invalid_argument("Error: Vectors must have the same size for subtraction");
    }

    std::vector<double> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); ++i) {
        result[i] = vector1[i] - vector2[i];
    }
    return result;
}

double vecNorm2(const std::vector<double>& vector) {
    double sum = 0.0;
    for (const auto& value : vector) {
        sum += value * value;
    }
    return sqrt(sum);
}

std::vector<double> getColumn(const std::vector<std::vector<double>>& matrix, size_t colIndex) {
    if (colIndex >= matrix[0].size() - 1) {
        throw std::out_of_range("Error: Invalid column index");
    }
    std::vector<double> column;
    for (const auto& row : matrix) {
        column.push_back(row[colIndex]);
    }
    return column;
}

size_t MaxAbs(const std::vector<double>& x, size_t i_init, size_t i_max, size_t i_step){
    if (i_init >= x.size()) {
        throw std::out_of_range("Error: Invalid index");
    }
    double max_v = std::abs(x[i_init]);
    size_t max_i = i_init;
    double v;
    for (size_t i=i_init+i_step; i < std::min(i_max, x.size()); i += i_step){
        v = std::abs(x[i]);
        if (v > max_v){
            max_v = v;
            max_i = i;
        }
    }
    return max_i;
}

size_t MaxAbs(const std::vector<double>& x, size_t i_init, size_t i_max){
    return MaxAbs(x, i_init, i_max, 1);
}
size_t MaxAbs(const std::vector<double>& x){
    return MaxAbs(x, 0, x.size());
}

std::vector<size_t> MaxAbsLine(const std::vector<std::vector<double>>& x, bool cumulative, size_t col_max){
    std::vector<size_t> i_list;
    size_t n = x.size();
    size_t col_init = 0;
    for (size_t i=0; i<n; i++){
        if (cumulative){
            col_init = i;
        }
        i_list.push_back(MaxAbs(x[i], col_init, col_max));
    }
    return i_list;
}
std::vector<size_t> MaxAbsLine(const std::vector<std::vector<double>>& x){
    return MaxAbsLine(x, false, x.size());
}

std::vector<size_t> MaxAbsCol(const std::vector<std::vector<double>>& x, bool cumulative, size_t col_max){
    std::vector<size_t> i_list;
    size_t i_init = 0;
    for (size_t i=0; i < std::min(col_max, x[0].size()); i++){
        if (cumulative){
            i_init = i;
        }
        i_list.push_back(MaxAbs(getColumn(x, i), i_init, x.size()));
    }
    return i_list;
}
std::vector<size_t> MaxAbsCol(const std::vector<std::vector<double>>& x){
    return MaxAbsCol(x, false, x[0].size());
}

std::vector<double> SolveGauss(
    const std::vector<std::vector<double>>& a_mat,
    const std::vector<double>& b_vec,
    bool with_scaling){

    if (a_mat.size() != b_vec.size() || a_mat.size() != a_mat[0].size()) {
        throw std::invalid_argument("Error: Invalid matrix or vector dimensions");
    }

    size_t n = a_mat.size();
    std::vector<size_t> n_line(n);
    std::vector<double> s;
    if (with_scaling){
        n_line = MaxAbsLine(a_mat);
        for (size_t i=0; i<n; i++){
            if (std::abs(a_mat[i][n_line[i]]) < eps){
                std::cerr << "Input matrix has a line of zeros. Response is non-unique! Cannot continue." << std::endl;
            }
            s.push_back(std::abs(a_mat[i][n_line[i]]));
        }
    }
    for (size_t i=0; i<n; i++){
        n_line[i] = i;
    }

    std::vector<std::vector<double>> a = combineMatrixAndVector(a_mat, b_vec);
    size_t p;
    double max_v;
    double v;
    double m;
    for (size_t i=0; i<n-1; i++){
        max_v = std::abs(a[n_line[i]][i]);
        if (with_scaling){
            max_v /= s[n_line[i]];
        }
        p = i;
        for (size_t j=i+1; j<n; j++){
            v = std::abs(a[n_line[j]][i]);
            if (with_scaling){
                v /= s[n_line[j]];;
            }
            if (v > max_v){
                max_v = v;
                p = j;
            }
        }
        if (std::abs(a[n_line[p]][i]) < eps){
            std::cerr << "Input matrix has a column of zeros. Response is non-unique! Cannot continue." << std::endl;
        } else if (n_line[i] != n_line[p]){
            std::swap(n_line[i], n_line[p]);
        }
        for (size_t j=i+1; j<n; j++){
            m = a[n_line[j]][i] / a[n_line[i]][i];
            for (size_t k=0; k<=n; k++){
                a[n_line[j]][k] = a[n_line[j]][k] - m * a[n_line[i]][k];
            }
        }
    }
    if (std::abs(a[n_line[n-1]][n-1]) < eps){
        std::cerr << "Equations are not lineally independent. Response is non-unique! Cannot continue." << std::endl;
    }

    std::vector<double> x(n);
    double s_xa;
    x[n-1] = a[n_line[n-1]][n] / a[n_line[n-1]][n-1];
    for (size_t i=n-1; i>0; i--){
        s_xa = 0.;
        for (size_t j=i+1; j<=n; j++){
            s_xa += a[n_line[i-1]][j-1] * x[j-1];
        }
        x[i-1] = (a[n_line[i-1]][n] - s_xa)/ a[n_line[i-1]][i-1];
    }
    return x;
}

std::vector<double> SolveGauss(
    const std::vector<std::vector<double>>& a_mat,
    const std::vector<double>& b_vec){
    return SolveGauss(a_mat, b_vec, true);
}

void check_result(
    const std::vector<std::vector<double>>& a_mat,
    const std::vector<double>& b_vec,
    std::vector<double>& x_out,
    const std::vector<double>& x_true,
    std::ostream& output){

    output << "Response" << std::endl;
    std::vector<double> error = vecDif(b_vec, matMult(a_mat, x_out));
    // printVector(x_out, output);
    // printVector(error, output);
    output << "  Maximum residual: " << std::abs(error[MaxAbs(error)]) << std::endl;
    output << "  ||residual||: " << vecNorm2(error) << std::endl;

    output << "Provided True Response" << std::endl;
    std::vector<double> error_true = vecDif(b_vec, matMult(a_mat, x_true));
    // printVector(x_true, output);
    // printVector(error_true, output);
    output << "  Maximum residual: " << std::abs(error_true[MaxAbs(error_true)]) << std::endl;
    output << "  ||residual||: " << vecNorm2(error_true) << std::endl;

    output << "Delta Responses" << std::endl;
    std::vector<double> x_delta = vecDif(x_true, x_out);
    // printVector(x_delta, output);
    output << "  ||delta||: " << vecNorm2(x_delta) << std::endl;
    output << "  Maximum difference: " << std::abs(x_delta[MaxAbs(x_delta)]) << std::endl;
    output << "  Maximum difference odd: " << std::abs(x_delta[MaxAbs(x_delta, 0, x_delta.size(), 2)]) << std::endl;
    output << "  Maximum difference even: " << std::abs(x_delta[MaxAbs(x_delta, 1, x_delta.size(), 2)]) << std::endl;
}

void test01(){
    std::vector<std::vector<double>> matrix = {
        {30., 591400.},
        {5.291, -6.130}
    };

    std::vector<double> vector = {591700., 46.78}; // Example vector
    std::vector<double> x_true = {10., 1.}; // Example vector

    std::vector<double> x_out;
    x_out = SolveGauss(matrix, vector, false);

    std::cout << "x = ";
    std::string sep = "{";
    for (const auto& v : x_out) {
        std::cout << sep<< v ;
        sep = ", ";
    }
    std::cout << "}" << std::endl;

    check_result(matrix, vector, x_out, x_true, std::cout);
}

void tests(std::ostream& output, bool is_2D, bool both_versions){
    std::string n;
    std::string n2 = "1";
    fs::path k_filename;
    fs::path f_filename;
    fs::path x_true_filename;
    fs::path x_est_filename;

    fs::path currentPath = fs::current_path();

    std::vector<std::vector<double>> k_matrix;
    std::vector<double> f_vector;
    std::vector<double> x_true;
    std::vector<double> x_est;
    std::vector<double> x_out;

    std::vector<int> a = {3, 5, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100};
    for (int i : a) {
        n = std::to_string(i);
        if (is_2D){
            n2 = n;
        }
        output << "Number of cells: " << n << " x " << n2 << std::endl;

        k_filename = currentPath / ("no_sync/matrix/k_" + n + "_" + n2 + ".csv");
        f_filename = currentPath / ("no_sync/matrix/f_" + n + "_" + n2 + ".csv");
        x_true_filename = currentPath / ("no_sync/matrix/x_" + n + "_" + n2 + ".csv");
        x_est_filename = currentPath / ("no_sync/matrix/x_prev_" + n + "_" + n2 + ".csv");

        std::cout << "Reading csv files (" << n << " x " << n2 << ")..";
        k_matrix = readCSV(k_filename.string());
        f_vector = flattenMatrix(readCSV(f_filename.string()));
        x_true = flattenMatrix(readCSV(x_true_filename.string()));
        x_est = flattenMatrix(readCSV(x_est_filename.string()));
        std::cout << " Done!" << std::endl;

        output << "  K matrix = " << k_matrix.size() << " x " << k_matrix[0].size() << std::endl;
        // output << "  F vector = " << f_vector.size() << " x 1" << std::endl;
        // output << "  X vector = " << x_true.size() << " x 1" << std::endl;
        // output << "  X est vector = " << x_est.size() << " x 1" << std::endl;

        output << "  Gauss Elimination with Partial Pivot and Scaling" << std::endl;
        auto start = high_resolution_clock::now();
        x_out = SolveGauss(k_matrix, f_vector, true);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        output << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;

        check_result(k_matrix, f_vector, x_out, x_true, output);
        output << std::endl;

        if (both_versions){
            output << "  Gauss Elimination with Partial Pivot" << std::endl;
            start = high_resolution_clock::now();
            x_out = SolveGauss(k_matrix, f_vector, false);
            stop = high_resolution_clock::now();
            duration = duration_cast<milliseconds>(stop - start);
            output << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;

            check_result(k_matrix, f_vector, x_out, x_true, output);
            output << std::endl;
        }

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

    // test01();
    std::ofstream outputFile("results_1D.txt");
    if (outputFile.is_open()) {
        tests(outputFile, false, false);  // 1D
        // tests(outputFile, true, false);  // 2D
        outputFile.close();
    } else {
        std::cerr << "Unable to open the file for writing.\n";
    }
}