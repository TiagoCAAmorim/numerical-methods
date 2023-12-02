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

class Info{
    public:
        Info();
        void add_value_int(std::string name, int value);
        void add_value_double(std::string name, double value);
        void reset_int();
        void reset_double();
        int get_value_int(std::string name);
        double get_value_double(std::string name);
    private:
        size_t position_i(std::string name);
        size_t position_d(std::string name);

        std::vector<std::string> names_i;
        std::vector<std::string> names_d;
        std::vector<int> values_i;
        std::vector<double> values_d;
};

Info::Info():
    names_i(),
    names_d(),
    values_i(),
    values_d()
    {};

size_t Info::position_i(std::string name){
    for (size_t j=0; j<names_i.size(); j++){
        if (names_i[j] == name){
            return j;
        }
    }
    return names_i.size();
}
size_t Info::position_d(std::string name){
    for (size_t j=0; j<names_d.size(); j++){
        if (names_d[j] == name){
            return j;
        }
    }
    return names_d.size();
}

void Info::add_value_int(std::string name, int value){
    size_t i = position_i(name);
    if (i < names_i.size()){
        names_i[i] = name;
        values_i[i] = value;
    } else{
        names_i.push_back(name);
        values_i.push_back(value);
    }
}
void Info::add_value_double(std::string name, double value){
    size_t i = position_d(name);
    if (i < names_d.size()){
        names_d[i] = name;
        values_d[i] = value;
    } else{
        names_d.push_back(name);
        values_d.push_back(value);
    }
}

void Info::reset_int(){
    names_i.clear();
    values_i.clear();
}
void Info::reset_double(){
    names_d.clear();
    values_d.clear();
}

int Info::get_value_int(std::string name){
    size_t i = position_i(name);
    if (i < names_i.size()){
        return values_i[i];
    }
    return -99999;
}
double Info::get_value_double(std::string name){
    size_t i = position_d(name);
    if (i < names_d.size()){
        return values_d[i];
    }
    return -99999.;
}

Info info;

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
double MaxAbsValue(const std::vector<double>& x){
    return std::abs(x[MaxAbs(x, 0, x.size())]);
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

std::vector<double> SolveSRS(
    const std::vector<std::vector<double>>& a_mat,
    const std::vector<double>& b_vec,
    const std::vector<double>& x0_vec,
    const double conv_tol,
    const size_t max_iterations,
    const double param_w,
    const bool is_jacobi){
//532
    std::vector<double> x_old = x0_vec;
    std::vector<double> x_new(a_mat.size());
    size_t n = a_mat.size();
    double a_x;
    double max_x;
    double max_dx;
    double conv;

    if (param_w <= 0.){
        std::cerr << "Parameter w must be greater than zero. Cannot continue" << std::endl;
    }

    for (size_t i=0; i < n; i++){
        if (std::abs(a_mat[i][i]) < eps){
            std::cerr << "Found zero in main diagonal. Reorganize data first. Cannot continue" << std::endl;
        }
    }

    for (size_t k=0; k < max_iterations; k++){
        for (size_t i=0; i < n; i++){
            a_x = 0;
            for (size_t j=0; j < n; j++){
                if (j==i){
                    // skip current variable
                } else if (is_jacobi || j > i){
                    a_x += a_mat[i][j] * x_old[j];
                } else{
                    a_x += a_mat[i][j] * x_new[j];
                }
            }
            x_new[i] = (1 - param_w) * x_old[i] + param_w / a_mat[i][i] * (b_vec[i] - a_x);
        }

        max_dx = MaxAbsValue(vecDif(x_new, x_old));
        max_x = MaxAbsValue(x_new);
        if (max_x < eps){
            max_x = MaxAbsValue(x_old);
            if (max_x < eps){
                conv = 0.;
                break;
            }
        }
        conv = max_dx / max_x;
        if (conv < conv_tol){
            break;
        }

        x_old = x_new;
        info.add_value_int("Iterations",static_cast<int>(k)+1);
    }

    info.add_value_double("Convergence",conv);
    return x_new;
}

std::vector<double> SolveJacobi(
    const std::vector<std::vector<double>>& a_mat,
    const std::vector<double>& b_vec,
    const std::vector<double>& x0_vec,
    const double conv_tol,
    const size_t max_iterations){
//517
    return SolveSRS(a_mat, b_vec, x0_vec, conv_tol, max_iterations, 1., true);
}

std::vector<double> SolveGaussSeidel(
    const std::vector<std::vector<double>>& a_mat,
    const std::vector<double>& b_vec,
    const std::vector<double>& x0_vec,
    const double conv_tol,
    const size_t max_iterations){
//520
    return SolveSRS(a_mat, b_vec, x0_vec, conv_tol, max_iterations, 1., false);
}

std::vector<double> SolveSRS(
    const std::vector<std::vector<double>>& a_mat,
    const std::vector<double>& b_vec,
    const std::vector<double>& x0_vec,
    const double conv_tol,
    const size_t max_iterations,
    const double param_w){
//520
    return SolveSRS(a_mat, b_vec, x0_vec, conv_tol, max_iterations, param_w, false);
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
        {10., -1., 2., 0.},
        {-1., 11., -1., 3.},
        {2., -1., 10., -1.},
        {0., 3., -1., 8.},
    };

    std::vector<double> vector = {6., 25., -11., 15.};
    std::vector<double> x_true = {1., 2., -1., 1.};
    std::vector<double> x_est  = {0., 0., 0., 0.};

    std::vector<double> x_out;

    std::cout << "Gauss" << std::endl;
    x_out = SolveGauss(matrix, vector, false);
    printVector(x_out, std::cout);
    check_result(matrix, vector, x_out, x_true, std::cout);

    std::cout << "Jacobi" << std::endl;
    x_out = SolveJacobi(matrix, vector, x_est, 1E-3, 10);
    printVector(x_out, std::cout);
    std::cout << "Iterations: " << info.get_value_int("Iterations") << std::endl;
    std::cout << "Convergence: " << info.get_value_double("Convergence") << std::endl;
    check_result(matrix, vector, x_out, x_true, std::cout);

    std::cout << "Gauss-Seidel" << std::endl;
    x_out = SolveGaussSeidel(matrix, vector, x_est, 1E-3, 10);
    printVector(x_out, std::cout);
    std::cout << "Iterations: " << info.get_value_int("Iterations") << std::endl;
    std::cout << "Convergence: " << info.get_value_double("Convergence") << std::endl;
    check_result(matrix, vector, x_out, x_true, std::cout);

    std::cout << "SRS - w = 0.5" << std::endl;
    x_out = SolveSRS(matrix, vector, x_est, 1E-3, 10, 0.5);
    printVector(x_out, std::cout);
    std::cout << "Iterations: " << info.get_value_int("Iterations") << std::endl;
    std::cout << "Convergence: " << info.get_value_double("Convergence") << std::endl;
    check_result(matrix, vector, x_out, x_true, std::cout);

    std::cout << "SRS - w = 1.5" << std::endl;
    x_out = SolveSRS(matrix, vector, x_est, 1E-3, 10, 1.5);
    printVector(x_out, std::cout);
    std::cout << "Iterations: " << info.get_value_int("Iterations") << std::endl;
    std::cout << "Convergence: " << info.get_value_double("Convergence") << std::endl;
    check_result(matrix, vector, x_out, x_true, std::cout);
}

void tests(
    std::ostream& output,
    const bool is_2D,
    const bool use_gauss,
    const double conv_tol,
    const size_t max_iterations,
    const double param_w1,
    const double param_w2){

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
        printVector(x_true, output);

        output << "  Jacobi" << std::endl;
        auto start = high_resolution_clock::now();
        x_out = SolveJacobi(k_matrix, f_vector, x_est, conv_tol, max_iterations);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        printVector(x_out, output);
        output << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;
        output << "Iterations: " << info.get_value_int("Iterations") << std::endl;
        output << "Convergence: " << info.get_value_double("Convergence") << std::endl;
        check_result(k_matrix, f_vector, x_out, x_true, output);
        output << std::endl;

        output << "  Gauss-Siedel" << std::endl;
        start = high_resolution_clock::now();
        x_out = SolveGaussSeidel(k_matrix, f_vector, x_est, conv_tol, max_iterations);
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        printVector(x_out, output);
        output << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;
        output << "Iterations: " << info.get_value_int("Iterations") << std::endl;
        output << "Convergence: " << info.get_value_double("Convergence") << std::endl;
        check_result(k_matrix, f_vector, x_out, x_true, output);
        output << std::endl;

        output << "  SRS - w = " << param_w1 << std::endl;
        start = high_resolution_clock::now();
        x_out = SolveSRS(k_matrix, f_vector, x_est, conv_tol, max_iterations, param_w1);
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        printVector(x_out, output);
        output << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;
        output << "Iterations: " << info.get_value_int("Iterations") << std::endl;
        output << "Convergence: " << info.get_value_double("Convergence") << std::endl;
        check_result(k_matrix, f_vector, x_out, x_true, output);
        output << std::endl;

        output << "  SRS - w = " << param_w2 << std::endl;
        start = high_resolution_clock::now();
        x_out = SolveSRS(k_matrix, f_vector, x_est, conv_tol, max_iterations, param_w2);
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        printVector(x_out, output);
        output << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;
        output << "Iterations: " << info.get_value_int("Iterations") << std::endl;
        output << "Convergence: " << info.get_value_double("Convergence") << std::endl;
        check_result(k_matrix, f_vector, x_out, x_true, output);
        output << std::endl;

        if (use_gauss){
            output << "  Gauss Elimination with Partial Pivot and Scaling" << std::endl;
            start = high_resolution_clock::now();
            x_out = SolveGauss(k_matrix, f_vector, true);
            stop = high_resolution_clock::now();
            duration = duration_cast<milliseconds>(stop - start);
            printVector(x_out, output);
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
    std::ofstream outputFile("results.txt");
    if (outputFile.is_open()) {
        tests(outputFile, true, false, 1E-5, 1000, 0.8, 1.2);
        outputFile.close();
    } else {
        std::cerr << "Unable to open the file for writing.\n";
    }
}