#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <iomanip>
// using namespace std;
// using namespace std::chrono;
// namespace fs = std::filesystem;

const double eps = 1e-17;

std::string d2str(double value, int precision){
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}
std::string d2str(double value){
    return d2str(value, 5);
}

std::string concatenate() {
    return "";
}
template <typename T, typename... Args>
std::string concatenate(T first, Args... args) {
    std::ostringstream oss;
    oss << first << concatenate(args...);
    return oss.str();
}

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

void print_matrix(std::vector<std::vector<double>> matrix, std::ostream& output){
    for (size_t r=0; r<matrix.size(); r++){
        std::string sep = "| ";
        for (size_t c=0; c<matrix[r].size(); c++){
            output << sep << matrix[r][c];
            sep = ", ";
        }
        output << " |" << std::endl;
    }
}

void print_vector(const std::string name, const std::vector<double>& vector, std::ostream& output, size_t first, size_t step){
    output << name;
    std::string sep = "{";
    for (size_t i=first; i<vector.size(); i+=step){
        output << sep << vector[i];
        sep = ", ";
    }
    output << "}";
}
void print_vector(const std::string name, const std::vector<double>& vector, std::ostream& output){
    print_vector(name, vector, output, 0, 1);
    output << std::endl;
}
void print_vector(const std::vector<double>& vector, std::ostream& output){
    print_vector("", vector, output);
    output << std::endl;
}

void initial_vector(std::vector<double>& vector, size_t n_rows, double init_value){
    vector.clear();
    vector.resize(n_rows, init_value);
}
void initial_vector(std::vector<double>& vector, size_t n_rows){
    initial_vector(vector, n_rows, 0.);
}

void initial_matrix(std::vector<std::vector<double>>& matrix, size_t n_rows, size_t n_cols, double init_value){
    matrix.resize(n_rows);

    for (auto& row : matrix) {
        row.resize(n_cols);
        std::fill(row.begin(), row.end(), init_value);
    }
}
void initial_matrix(std::vector<std::vector<double>>& matrix, size_t n_rows, size_t n_cols){
    initial_matrix(matrix, n_cols, n_rows, 0.);
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

std::vector<double> vecAdd(
    const std::vector<double>& vector1,
    const std::vector<double>& vector2) {

    if (vector1.size() != vector2.size()) {
        throw std::invalid_argument("Error: Vectors must have the same size for subtraction");
    }

    std::vector<double> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); ++i) {
        result[i] = vector1[i] + vector2[i];
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
            if (m != 0.){
                for (size_t k=0; k<=n; k++){
                    a[n_line[j]][k] = a[n_line[j]][k] - m * a[n_line[i]][k];
                }
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
    info.add_value_int("GaussCalls",info.get_value_int("GaussCalls")+1);
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


class RelPerm{
    public:
        RelPerm();

        void set_swi(double value);
        void set_swc(double value);
        void set_sor(double value);

        void set_nw(double value);
        void set_no(double value);
        void set_krw_endp(double value);
        void set_kro_endp(double value);

        double get_swi();
        double get_krw(double sw);
        double get_kro(double sw);
        double get_dkrw(double sw);
        double get_dkro(double sw);
    private:
        double get_swd(double sw);
        double get_sod(double sw);

        double swi;
        double swc;
        double sor;

        double nw;
        double no;
        double krw_endp;
        double kro_endp;

};
RelPerm::RelPerm():
    swi(0.),
    swc(0.),
    sor(0.),
    nw(1.),
    no(1.),
    krw_endp(1.),
    kro_endp(1.)
    {};

void RelPerm::set_swi(double value){
    swi = value;
}
void RelPerm::set_swc(double value){
    swc = value;
}
void RelPerm::set_sor(double value){
    sor = value;
}

void RelPerm::set_nw(double value){
    nw = value;
}
void RelPerm::set_no(double value){
    no = value;
}
void RelPerm::set_krw_endp(double value){
    krw_endp = value;
}
void RelPerm::set_kro_endp(double value){
    kro_endp = value;
}

double RelPerm::get_swi(){
    return swi;
}
double RelPerm::get_swd(double sw){
    if (sw <= swc){
        return 0.;
    } else if (sw >= (1. - sor)){
        return 1.;
    } else{
        return (sw - swc) / (1. - sor - swc);
    }
}
double RelPerm::get_sod(double sw){
    if (sw <= swi){
        return 1.;
    } else if (sw >= (1. - sor)){
        return 0.;
    } else{
        return 1. - (sw - swi) / (1. - sor - swi);
    }
}

double RelPerm::get_krw(double sw){
    if (sw <= swc){
        return 0.;
    } else if (sw >= (1. - sor)){
        return krw_endp + (1. - krw_endp) * (sw - (1. - sor)) / sor;
    } else{
        return krw_endp * pow(get_swd(sw), nw);
    }
}
double RelPerm::get_kro(double sw){
    if (sw <= swi){
        return kro_endp;
    } else if (sw >= (1. - sor)){
        return 0.;
    } else{
        return kro_endp * pow(get_sod(sw), no);
    }
}

double RelPerm::get_dkrw(double sw){
    if (sw <= swc){
        return 0.;
    } else if (sw >= (1. - sor)){
        return (1. - krw_endp) / sor;
    } else{
        return krw_endp * nw * pow(get_swd(sw), nw - 1.) / (1. - sor - swc);
    }
}
double RelPerm::get_dkro(double sw){
    if (sw <= swi){
        return 0.;
    } else if (sw >= (1. - sor)){
        return 0.;
    } else{
        return (-1.) * kro_endp * no * pow(get_sod(sw), no - 1.) / (1. - sor - swi);
    }
}


class SimModel{
    public:
        SimModel();

        // ResultProcessor(std::ostream& output = std::cout) : outputStream(output) {}

        void set_output(std::ostream& value);

        void set_t_end(double value);
        void set_max_dsw(double value);
        void set_max_dpr(double value);
        void set_max_dt(double value);
        void set_min_dt(double value);
        void set_conv_tol(double value);
        void set_max_iter(size_t value);
        void set_use_nr(bool value);

        void set_ni(size_t value);
        void set_nj(size_t value);

        void set_dh(double value);
        void set_dz(double value);
        void set_phi(double value);
        void set_perm(double value);
        void set_p_init(double value);

        void set_bo(double value);
        void set_bw(double value);
        void set_uo(double value);
        void set_uw(double value);

        void set_qwi(double value);
        void set_pwf(double value);
        void set_rw(double value);
        void set_skin(double value);

        void run_simulation(double dt);

        std::vector<double> get_t_vec();
        std::vector<double> get_dt_vec();
        std::vector<double> get_max_dsw_vec();
        std::vector<double> get_max_dpr_vec();
        std::vector<double> get_qo_vec();
        std::vector<double> get_qw_vec();
        std::vector<double> get_pinj_vec();

        RelPerm kr;
    private:
        size_t get_cell_n(size_t i, size_t j);
        double get_cell_pr(std::vector<double> x, size_t c);
        double get_cell_pr(std::vector<double> x, size_t i, size_t j);
        double get_cell_sw(std::vector<double> x, size_t c);
        double get_cell_sw(std::vector<double> x, size_t i, size_t j);

        double get_upstream_sw(std::vector<double> x, size_t c1, size_t c2);
        double get_upstream_sw(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2);

        double get_kro(std::vector<double> x, size_t c);
        double get_kro(std::vector<double> x, size_t c1, size_t c2);
        double get_kro(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2);
        double get_krw(std::vector<double> x, size_t c);
        double get_krw(std::vector<double> x, size_t c1, size_t c2);
        double get_krw(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2);

        double get_dkro(std::vector<double> x, size_t c1, size_t c2);
        double get_dkro(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2);
        double get_dkrw(std::vector<double> x, size_t c1, size_t c2);
        double get_dkrw(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2);

        double get_tr(size_t c1, size_t c2);
        double get_tr(size_t i1, size_t j1, size_t i2, size_t j2);
        double get_tro(std::vector<double> x, size_t c1, size_t c2);
        double get_tro(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2);
        double get_trw(std::vector<double> x, size_t c1, size_t c2);
        double get_trw(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2);

        double get_dtro(std::vector<double> x, size_t c1, size_t c2);
        double get_dtro(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2);
        double get_dtrw(std::vector<double> x, size_t c1, size_t c2);
        double get_dtrw(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2);

        double get_dtrodP(std::vector<double> x, size_t c1, size_t c2);
        double get_dtrwdP(std::vector<double> x, size_t c1, size_t c2);

        void build_K(std::vector<double> x);

        void build_F(std::vector<double> x);
        void build_Jac(std::vector<double> x);

        void initialize();
        double get_convergence(std::vector<double>& x_old, std::vector<double>& x_new);
        bool solve_fixed();
        bool solve_newton_raphson();

        double get_max_dpr();
        double get_max_dsw();
        bool advance_one_time_step();
        void save_result();

        void log_msg(std::string message);

        std::ostream& output;

        double t_end;
        double max_dsw;
        double max_dpr;
        double max_dt;
        double min_dt;
        double conv_tol;
        size_t max_iter;
        bool use_nr;
        double dt_curr;

        size_t ni;
        size_t nj;

        double dh;
        double di;
        double dj;
        double dz;
        double phi;
        double perm;
        double p_init;

        double bo;
        double bw;
        double uo;
        double uw;

        double qwi;
        double pwf;
        double rw;
        double skin;
        double wi;

        std::vector<std::vector<double>> k_mat;
        std::vector<double> f_vec;
        std::vector<std::vector<double>> jac;
        std::vector<double> x_curr;
        std::vector<double> x_next;
        std::vector<std::vector<double>> x_list;

        std::vector<double> t_vec;
        std::vector<double> dt_vec;
        std::vector<double> max_dsw_vec;
        std::vector<double> max_dpr_vec;

        std::vector<double> qo_vec;
        std::vector<double> qw_vec;
        std::vector<double> pinj_vec;
};

double unit_conv = 0.00852702; // units: bar, mD, cP, m, m3/d

SimModel::SimModel():
    kr(),

    output(std::cout),

    t_end(10.),
    max_dsw(0.005),
    max_dpr(1.),
    max_dt(10.),
    min_dt(0.01),
    conv_tol(1.E-5),
    max_iter(50),
    use_nr(false),
    dt_curr(0.1),

    ni(5),
    nj(1),
    dh(100.),
    di(100.),
    dj(100.),
    dz(100.),
    phi(0.01),
    perm(100.),
    p_init(250.),

    bo(1.2),
    bw(1.),
    uo(0.7),
    uw(1.),

    qwi(300.),
    pwf(200.),
    rw(0.1),
    skin(0.),
    wi(1.),

    k_mat(),
    f_vec(),
    jac(),
    x_curr(),
    x_next(),
    x_list(),

    t_vec(),
    dt_vec(),
    max_dsw_vec(),
    max_dpr_vec(),

    qo_vec(),
    qw_vec(),
    pinj_vec()
    {};

void SimModel::set_output(std::ostream& value){
    output.rdbuf(value.rdbuf());
}
void SimModel::log_msg(std::string message){
    output << message << std::endl;
}

void SimModel::set_t_end(double value){
    t_end = value;
}
void SimModel::set_max_dsw(double value){
    max_dsw = value;
}
void SimModel::set_max_dpr(double value){
    max_dpr = value;
}
void SimModel::set_max_dt(double value){
    max_dt = value;
}
void SimModel::set_min_dt(double value){
    min_dt = value;
}
void SimModel::set_conv_tol(double value){
    conv_tol = value;
}
void SimModel::set_max_iter(size_t value){
    max_iter = value;
}
void SimModel::set_use_nr(bool value){
    use_nr = value;
}
void SimModel::set_ni(size_t value){
    ni = value;
}
void SimModel::set_nj(size_t value){
    nj = value;
}
void SimModel::set_dh(double value){
    dh = value;
}
void SimModel::set_dz(double value){
    dz = value;
}
void SimModel::set_phi(double value){
    phi = value;
}
void SimModel::set_perm(double value){
    perm = value;
}
void SimModel::set_p_init(double value){
    p_init = value;
}
void SimModel::set_bo(double value){
    bo = value;
}
void SimModel::set_bw(double value){
    bw = value;
}
void SimModel::set_uo(double value){
    uo = value;
}
void SimModel::set_uw(double value){
    uw = value;
}
void SimModel::set_qwi(double value){
    qwi = value;
}
void SimModel::set_pwf(double value){
    pwf = value;
}
void SimModel::set_rw(double value){
    rw = value;
}
void SimModel::set_skin(double value){
    skin = value;
}

std::vector<double> SimModel::get_t_vec(){
    return t_vec;
}
std::vector<double> SimModel::get_dt_vec(){
    return dt_vec;
}
std::vector<double> SimModel::get_max_dsw_vec(){
    return max_dsw_vec;
}
std::vector<double> SimModel::get_max_dpr_vec(){
    return max_dpr_vec;
}
std::vector<double> SimModel::get_qo_vec(){
    return qo_vec;
}
std::vector<double> SimModel::get_qw_vec(){
    return qw_vec;
}
std::vector<double> SimModel::get_pinj_vec(){
    return pinj_vec;
}


size_t SimModel::get_cell_n(size_t i, size_t j){
    return i + 1 + j * ni;
}
double SimModel::get_cell_pr(std::vector<double> x, size_t c){
    return x[2 * c - 2];
}
double SimModel::get_cell_pr(std::vector<double> x, size_t i, size_t j){
    size_t c = get_cell_n(i,j);
    return get_cell_pr(x, c);
}
double SimModel::get_cell_sw(std::vector<double> x, size_t c){
    return x[2 * c - 1];
}
double SimModel::get_cell_sw(std::vector<double> x, size_t i, size_t j){
    size_t c = get_cell_n(i,j);
    return get_cell_sw(x, c);
}

double SimModel::get_upstream_sw(std::vector<double> x, size_t c1, size_t c2){
    if (c1==c2){
        return get_cell_sw(x,c1);
    } else if (get_cell_pr(x,c1) > get_cell_pr(x,c2)){
        return get_cell_sw(x,c1);
    } else{
        return get_cell_sw(x,c2);
    }
}
double SimModel::get_upstream_sw(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2){
    size_t c1 = get_cell_n(i1,j1);
    size_t c2 = get_cell_n(i2,j2);
    return get_upstream_sw(x, c1, c2);
}

double SimModel::get_kro(std::vector<double> x, size_t c){
    return get_kro(x, c, c);
}
double SimModel::get_kro(std::vector<double> x, size_t c1, size_t c2){
    double swi = get_upstream_sw(x,c1,c2);
    return kr.get_kro(swi);
}
double SimModel::get_kro(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2){
    size_t c1 = get_cell_n(i1,j1);
    size_t c2 = get_cell_n(i2,j2);
    return get_kro(x, c1, c2);
}
double SimModel::get_krw(std::vector<double> x, size_t c){
    return get_krw(x, c, c);
}
double SimModel::get_krw(std::vector<double> x, size_t c1, size_t c2){
    double swi = get_upstream_sw(x,c1,c2);
    return kr.get_krw(swi);
}
double SimModel::get_krw(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2){
    size_t c1 = get_cell_n(i1,j1);
    size_t c2 = get_cell_n(i2,j2);
    return get_krw(x, c1, c2);
}
double SimModel::get_dkro(std::vector<double> x, size_t c1, size_t c2){
    double swi = get_upstream_sw(x,c1,c2);
    return kr.get_dkro(swi);
}
double SimModel::get_dkro(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2){
    size_t c1 = get_cell_n(i1,j1);
    size_t c2 = get_cell_n(i2,j2);
    return get_dkro(x, c1, c2);
}
double SimModel::get_dkrw(std::vector<double> x, size_t c1, size_t c2){
    double swi = get_upstream_sw(x,c1,c2);
    return kr.get_dkrw(swi);
}
double SimModel::get_dkrw(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2){
    size_t c1 = get_cell_n(i1,j1);
    size_t c2 = get_cell_n(i2,j2);
    return get_dkrw(x, c1, c2);
}

double SimModel::get_tr(size_t c1, size_t c2){
    if (std::max(c1,c2) - std::min(c1,c2) == ni){
        return perm * di * dz / dj;
    }
    if (std::max(c1,c2) - std::min(c1,c2) == 1){
        return perm * dj * dz / di;
    }
    return 0.;
}
double SimModel::get_tr(size_t i1, size_t j1, size_t i2, size_t j2){
    size_t c1 = get_cell_n(i1,j1);
    size_t c2 = get_cell_n(i2,j2);
    return get_tr(c1, c2);
}

double SimModel::get_tro(std::vector<double> x, size_t c1, size_t c2){
    double tr = get_tr(c1,c2);
    if (tr != 0.){
        tr *= unit_conv * get_kro(x,c1,c2) / (bo * uo);
    }
    return tr;
}
double SimModel::get_tro(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2){
    size_t c1 = get_cell_n(i1,j1);
    size_t c2 = get_cell_n(i2,j2);
    return get_tro(x, c1, c2);
}
double SimModel::get_trw(std::vector<double> x, size_t c1, size_t c2){
    double tr = get_tr(c1,c2);
    if (tr != 0.){
        tr *= unit_conv * get_krw(x,c1,c2) / (bw * uw);
    }
    return tr;
}
double SimModel::get_trw(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2){
    size_t c1 = get_cell_n(i1,j1);
    size_t c2 = get_cell_n(i2,j2);
    return get_trw(x, c1, c2);
}

double SimModel::get_dtro(std::vector<double> x, size_t c1, size_t c2){
    double dtr = get_tr(c1,c2);
    if (dtr != 0.){
        dtr *= unit_conv * get_dkro(x,c1,c2) / (bo * uo);
    }
    return dtr;
}
double SimModel::get_dtro(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2){
    size_t c1 = get_cell_n(i1,j1);
    size_t c2 = get_cell_n(i2,j2);
    return get_dtro(x, c1, c2);
}
double SimModel::get_dtrw(std::vector<double> x, size_t c1, size_t c2){
    double dtr = get_tr(c1,c2);
    if (dtr != 0.){
        dtr *= unit_conv * get_dkrw(x,c1,c2) / (bw * uw);
    }
    return dtr;
}
double SimModel::get_dtrw(std::vector<double> x, size_t i1, size_t j1, size_t i2, size_t j2){
    size_t c1 = get_cell_n(i1,j1);
    size_t c2 = get_cell_n(i2,j2);
    return get_dtrw(x, c1, c2);
}

double SimModel::get_dtrodP(std::vector<double> x, size_t c1, size_t c2){
    double p1 = get_cell_pr(x, c1);
    double p2 = get_cell_pr(x, c2);
    if (p1 > p2){
        return get_dtro(x,c1,c2) * (p2-p1);
    } else{
        return 0.;
    }
}
double SimModel::get_dtrwdP(std::vector<double> x, size_t c1, size_t c2){
    double p1 = get_cell_pr(x, c1);
    double p2 = get_cell_pr(x, c2);
    if (p1 > p2){
        return get_dtrw(x,c1,c2) * (p2-p1);
    } else{
        return 0.;
    }
}

void SimModel::build_K(std::vector<double> x){
    initial_matrix(k_mat, 2*ni*nj, 2*ni*nj);
    std::vector<size_t> c_list;
    size_t nc;
    double tro;
    double trw;
    size_t c;

    double vp_dt = di * dj * dz * phi / dt_curr;
    double vp_dt_bo = vp_dt / bo;
    double vp_dt_bw = (-1.) * vp_dt / bw;
    for (size_t j=0; j<nj; j++){
        for (size_t i=0; i<ni; i++){
            c = get_cell_n(i,j);
            c_list.clear();
            if (i > 0){
                c_list.push_back(get_cell_n(i-1, j));
            }
            if (j > 0){
                c_list.push_back(get_cell_n(i, j-1));
            }
            if (i < (ni-1)){
                c_list.push_back(get_cell_n(i+1, j));
            }
            if (j < (nj-1)){
                c_list.push_back(get_cell_n(i, j+1));
            }
            nc = c_list.size();
            for (size_t k=0; k<nc; k++){
                tro = get_tro(x, c, c_list[k]);
                trw = get_trw(x, c, c_list[k]);
                k_mat[2*c-2][2*c-2] -= tro;
                k_mat[2*c-1][2*c-2] -= trw;
                k_mat[2*c-2][2*c_list[k]-2] = tro;
                k_mat[2*c-1][2*c_list[k]-2] = trw;
            }
            k_mat[2*c-2][2*c-1] = vp_dt_bo;
            k_mat[2*c-1][2*c-1] = vp_dt_bw;
        }
    }
    c = get_cell_n(ni-1,nj-1);
    k_mat[2*c - 2][2*c - 2] -= wi * get_kro(x, c) / (bo * uo);
    k_mat[2*c - 1][2*c - 2] -= wi * get_krw(x, c) / (bw * uw);
}

void SimModel::build_F(std::vector<double> x){
    initial_vector(f_vec, 2*ni*nj);
    size_t c;

    double vp_dt = di * dj * dz * phi / dt_curr;
    double vp_dt_bo = vp_dt / bo;
    double vp_dt_bw = (-1.) * vp_dt / bw;
    double sw;
    for (size_t j=0; j<nj; j++){
        for (size_t i=0; i<ni; i++){
            c = get_cell_n(i,j);
            sw = get_cell_sw(x, c);
            f_vec[2*c-2] = vp_dt_bo * sw;
            f_vec[2*c-1] = vp_dt_bw * sw;
        }
    }
    f_vec[1] -= qwi;
    c = get_cell_n(ni-1,nj-1);
    f_vec[2*c - 2] -= wi * get_kro(x, c) / (bo * uo) * pwf;
    f_vec[2*c - 1] -= wi * get_krw(x, c) / (bw * uw) * pwf;
}

void SimModel::build_Jac(std::vector<double> x){
    jac = k_mat;

    std::vector<size_t> c_list;
    size_t nc;
    double dtro;
    double dtrw;
    size_t c;

    for (size_t j=0; j<nj; j++){
        for (size_t i=0; i<ni; i++){
            c = get_cell_n(i,j);
            c_list.clear();
            if (i > 0){
                c_list.push_back(get_cell_n(i-1, j));
            }
            if (j > 0){
                c_list.push_back(get_cell_n(i, j-1));
            }
            if (i < (ni-1)){
                c_list.push_back(get_cell_n(i+1, j));
            }
            if (j < (nj-1)){
                c_list.push_back(get_cell_n(i, j+1));
            }
            nc = c_list.size();
            for (size_t k=0; k<nc; k++){
                dtro = get_dtrodP(x, c, c_list[k]);
                dtrw = get_dtrwdP(x, c, c_list[k]);
                jac[2*c-2][2*c-1] += dtro;
                jac[2*c-1][2*c-1] += dtrw;
                jac[2*c_list[k]-2][2*c-1] -= dtro;
                jac[2*c_list[k]-2][2*c-1] -= dtrw;
            }
        }
    }
    c = get_cell_n(ni-1,nj-1);
    double p_prod = get_cell_pr(x, ni-1, nj-1);
    jac[2*c - 2][2*c - 1] -= wi * get_dkro(x, c, c) / (bo * uo) * (p_prod - pwf);
    jac[2*c - 1][2*c - 1] -= wi * get_dkrw(x, c, c) / (bw * uw) * (p_prod - pwf);
}

void SimModel::initialize(){
    if (ni < 1){
        std::cerr << "ni must be a positive integer." << std::endl;
    }
    if (nj < 1){
        std::cerr << "nj must be a positive integer." << std::endl;
    }
    di = dh;
    dj = dh;

    double a = dj / di;
    double ro = di * exp((-1.)*(a * M_PI - log(a))/(1. + a * a));
    wi = unit_conv * 2. * M_PI * perm * dz / (log(ro/rw) + skin);

    t_vec.clear();
    t_vec.push_back(0.);
    x_list.clear();

    double swi = kr.get_swi();
    x_curr.clear();
    for (size_t j=0; j<nj; j++){
        for (size_t i=0; i<ni; i++){
            x_curr.push_back(p_init);
            x_curr.push_back(swi);
        }
    }

    x_next = x_curr;
    x_list.push_back(x_curr);

    dt_vec.clear();
    max_dsw_vec.clear();
    max_dpr_vec.clear();
    qo_vec.clear();
    qw_vec.clear();
    pinj_vec.clear();

    dt_vec.push_back(0.);
    max_dsw_vec.push_back(0.);
    max_dpr_vec.push_back(0.);
    qo_vec.push_back(0.);
    qw_vec.push_back(0.);
    pinj_vec.push_back(0.);
}

double SimModel::get_convergence(std::vector<double>& x_old, std::vector<double>& x_new){
    double max_conv = 0.;
    double conv;
    double delta;
    for (size_t i=0; i<x_old.size(); i++){
        delta = std::abs(x_old[i] - x_new[i]);
        if (std::abs(x_new[i]) > 1.E-15){
            conv = delta / std::abs(x_new[i]);
        } else if (std::abs(x_old[i]) > 1.E-15){
            conv = delta / std::abs(x_old[i]);
        } else{
            conv = 0.;
        }
        if (conv > max_conv){
            max_conv = conv;
        }
    }
    return max_conv;
}

bool SimModel::solve_fixed(){
    std::vector<double> x_old;
    std::vector<double> x_new = x_curr;
    size_t k=0;
    double conv = 1.e5;
    while ((conv > conv_tol) && (k < max_iter)){
        x_old = x_new;
        build_K(x_new);
        build_F(x_curr);
        x_new = SolveGauss(k_mat, f_vec);
        conv = get_convergence(x_old, x_new);
        log_msg(concatenate("   k=",k,", conv=",conv));
        k++;
    }
    x_next = x_new;
    return (conv <= conv_tol);
}
bool SimModel::solve_newton_raphson(){
    std::vector<double> x_old;
    std::vector<double> dx;
    std::vector<double> x_new = x_curr;
    std::vector<std::vector<double>> kx;
    std::vector<double> r;
    size_t k=0;
    double conv = 1.e5;
    while ((conv > conv_tol) && (k < max_iter)){
        x_old = x_new;
        build_K(x_new);
        build_F(x_curr);
        r = matMult(k_mat, x_new);
        // r = flattenMatrix(kx);
        r = vecDif(f_vec,r);
        build_Jac(x_new);
        dx = SolveGauss(jac, r);
        x_new = vecAdd(x_old, dx);
        conv = get_convergence(x_old, x_new);
        log_msg(concatenate("   k=",k,", conv=",conv));
        k++;
    }
    x_next = x_new;
    return (conv <= conv_tol);
}
double SimModel::get_max_dpr(){
    double dpr_max = 0.;
    double dpr;
    size_t c;
    for (size_t j=0; j<nj; j++){
        for (size_t i=0; i<ni; i++){
            c = get_cell_n(i,j);
            dpr = get_cell_pr(x_curr, c);
            dpr -= get_cell_pr(x_next, c);
            dpr = std::abs(dpr);
            if (dpr > dpr_max){
                dpr_max = dpr;
            }
        }
    }
    return dpr_max;
}
double SimModel::get_max_dsw(){
    double dsw_max = 0.;
    double dsw;
    size_t c;
    for (size_t j=0; j<nj; j++){
        for (size_t i=0; i<ni; i++){
            c = get_cell_n(i,j);
            dsw = get_cell_sw(x_curr, c);
            dsw -= get_cell_sw(x_next, c);
            dsw = std::abs(dsw);
            if (dsw > dsw_max){
                dsw_max = dsw;
            }
        }
    }
    return dsw_max;
}
bool SimModel::advance_one_time_step(){
    bool ok;
    if (use_nr){
        ok = solve_newton_raphson();
    } else{
        ok = solve_fixed();
    }
    double m_dpr = get_max_dpr();
    double m_dsw = get_max_dsw();
    // log_msg(concatenate("   MaxDpr=",m_dpr,", MaxDsw=",m_dsw));
    if (ok){
        if (m_dpr <= max_dpr){
            if (m_dsw <= max_dsw){
                return true;
            }
        }
    }
    return false;
}

void SimModel::run_simulation(double dt){
    initialize();
    double t = 0.;
    bool ok;
    dt_curr = dt;
    info.add_value_int("GaussCalls",0);
    while (t < t_end){
        // log_msg(concatenate("[",d2str(t,3),"], Time-step=",d2str(dt_curr,3)));
        ok = advance_one_time_step();
        if (ok || dt_curr <= min_dt){
            t += dt_curr;
            save_result();
            dt_curr *= 1.2;
            if (dt_curr > max_dt){
                dt_curr = max_dt;
            }
        } else{
            dt_curr = dt_curr / 2.;
            if (dt_curr < min_dt){
                dt_curr = min_dt;
            }
            // log_msg("["+d2str(t, 3) + "] Could not converge. New time-step: "+d2str(dt_curr,3)+".");
        }
        dt_curr = std::min(dt_curr, t_end - t);
    }
}

void SimModel::save_result(){
    t_vec.push_back(t_vec.back() + dt_curr);
    dt_vec.push_back(dt_curr);

    double m_dpr = get_max_dpr();
    double m_dsw = get_max_dsw();
    max_dpr_vec.push_back(m_dpr);
    max_dsw_vec.push_back(m_dsw);

    size_t c = get_cell_n(ni-1, nj-1);
    double qo = wi * get_kro(x_next, c) / (bo * uo) * (get_cell_pr(x_next, c) - pwf);
    double qw = wi * get_krw(x_next, c) / (bw * uw) * (get_cell_pr(x_next, c) - pwf);
    qo_vec.push_back(qo);
    qw_vec.push_back(qw);

    pinj_vec.push_back(get_cell_pr(x_next, 0, 0));

    x_list.push_back(x_next);
    x_curr = x_next;

    log_msg(concatenate("[",d2str(t_vec.back(),3),"], dt=",dt_curr,", MaxDpr=",m_dpr,", MaxDsw=",m_dsw,", Qo=",d2str(qo,2),", Qw=",d2str(qw,2),", Gauss=",info.get_value_int("GaussCalls")));
}


void test(
    std::ostream& output,
    const bool is_2D,
    const bool is_NR,
    const double max_dpr,
    const double max_dsw,
    const double max_dt,
    const double min_dt,
    const double conv_tol,
    const size_t max_iterations){

    SimModel model;

    model.set_output(output);
    model.set_t_end(5.*365.25);
    model.set_max_dsw(max_dsw);
    model.set_max_dpr(max_dpr);
    model.set_max_dt(max_dt);
    model.set_min_dt(min_dt);
    model.set_conv_tol(conv_tol);
    model.set_max_iter(max_iterations);
    model.set_use_nr(is_NR);

    size_t ni = 10;
    model.set_ni(ni);
    if (is_2D){
        model.set_nj(ni);
    } else{
        model.set_nj(1);
    }
    model.set_dh(50.);
    model.set_dz(30.);
    model.set_phi(0.15);
    model.set_perm(1000.);
    model.set_p_init(340.);
    model.set_bo(1.01);
    model.set_bw(1.0);
    model.set_uo(130.);
    model.set_uw(1.0);

    model.kr.set_swi(0.20);
    model.kr.set_swc(0.20);
    model.kr.set_sor(0.15);
    model.kr.set_nw(2.0);
    model.kr.set_no(3.0);
    model.kr.set_kro_endp(1.00);
    model.kr.set_krw_endp(0.60);

    model.set_qwi(350.);
    model.set_pwf(330.);
    model.set_rw(4.*2.54/100.);
    model.set_skin(0.0);

    model.run_simulation(0.1);
}

int main(){
    #ifdef _WIN32
        system("cls");
    #elif __unix__
        system("clear");
    #else
        std::cout << "Running on an unknown system" << std::endl;
    #endif

    // std::ofstream outputFile1("results_Fixed_1D.txt");
    // if (outputFile1.is_open()) {
    //     test(outputFile1, false, false, 1., 0.005, 10., 0.1, 1.E-5, 20);
    //     std::cout << "End of simulation" << std::endl;
    //     outputFile1.close();
    // } else {
    //     std::cerr << "Unable to open the file for writing.\n";
    // }

    // std::ofstream outputFile2("results_NR_1D.txt");
    // if (outputFile2.is_open()) {
    //     test(outputFile2, false, true, 1., 0.005, 10., 0.1, 1.E-5, 20);
    //     std::cout << "End of simulation" << std::endl;
    //     outputFile2.close();
    // } else {
    //     std::cerr << "Unable to open the file for writing.\n";
    // }

    std::ofstream outputFile3("results_Fixed_2D.txt");
    if (outputFile3.is_open()) {
        test(outputFile3, true, false, 1., 0.005, 10., 0.1, 1.E-5, 20);
        std::cout << "End of simulation" << std::endl;
        outputFile3.close();
    } else {
        std::cerr << "Unable to open the file for writing.\n";
    }

    std::ofstream outputFile4("results_NR_2D.txt");
    if (outputFile4.is_open()) {
        test(outputFile4, true, true, 1., 0.005, 10., 0.1, 1.E-5, 20);
        std::cout << "End of simulation" << std::endl;
        outputFile4.close();
    } else {
        std::cerr << "Unable to open the file for writing.\n";
    }
}