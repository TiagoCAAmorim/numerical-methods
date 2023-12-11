#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>
#include <chrono>
#include <cmath>
#include <algorithm>

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

void initial_vector(std::vector<double>& vector, int n_rows, double init_value){
    vector.clear();
    vector.resize(n_rows, init_value);
}
void initial_vector(std::vector<double>& vector, int n_rows){
    initial_vector(vector, n_rows, 0.);
}

void initial_matrix(std::vector<std::vector<double>>& matrix, int n_rows, int n_cols, double init_value){
    matrix.resize(n_rows);

    for (auto& row : matrix) {
        row.resize(n_cols);
        std::fill(row.begin(), row.end(), init_value);
    }
}
void initial_matrix(std::vector<std::vector<double>>& matrix, int n_rows, int n_cols){
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

        void set_t_end(double value);
        void set_max_dsw(double value);
        void set_max_dpr(double value);
        void set_max_dt(double value);
        void set_min_dt(double value);
        void set_conv_tol(double value);
        void set_use_nr(bool value);

        void set_ni(int value);
        void set_nj(int value);

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
    private:
        void initialize();

        int get_cell_n(int i, int j);
        double get_cell_pr(std::vector<double> x, int c);
        double get_cell_pr(std::vector<double> x, int i, int j);
        double get_cell_sw(std::vector<double> x, int c);
        double get_cell_sw(std::vector<double> x, int i, int j);

        double get_upstream_sw(std::vector<double> x, int c1, int c2);
        double get_upstream_sw(std::vector<double> x, int i1, int j1, int i2, int j2);
        double get_kro(std::vector<double> x, int c);
        double get_kro(std::vector<double> x, int c1, int c2);
        double get_kro(std::vector<double> x, int i1, int j1, int i2, int j2);
        double get_krw(std::vector<double> x, int c);
        double get_krw(std::vector<double> x, int c1, int c2);
        double get_krw(std::vector<double> x, int i1, int j1, int i2, int j2);
        double get_dkro(std::vector<double> x, int c1, int c2);
        double get_dkro(std::vector<double> x, int i1, int j1, int i2, int j2);
        double get_dkrw(std::vector<double> x, int c1, int c2);
        double get_dkrw(std::vector<double> x, int i1, int j1, int i2, int j2);
        double get_tr(int c1, int c2);
        double get_tr(int i1, int j1, int i2, int j2);
        double get_tro(std::vector<double> x, int c1, int c2);
        double get_tro(std::vector<double> x, int i1, int j1, int i2, int j2);
        double get_trw(std::vector<double> x, int c1, int c2);
        double get_trw(std::vector<double> x, int i1, int j1, int i2, int j2);

        void advance_one_time_step();
        void save_result();

        void build_K(std::vector<double> x, double dt);
        void build_K();

        void build_F(std::vector<double> x, double dt);
        void build_F();
        void build_Jac();

        double t_end;
        double max_dsw;
        double max_dpr;
        double max_dt;
        double min_dt;
        double conv_tol;
        bool use_nr;
        double dt_curr;

        int ni;
        int nj;

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
        RelPerm kr;

        double qwi;
        double pwf;
        double rw;
        double skin;
        double wi;

        std::vector<double> t_vec;
        std::vector<std::vector<double>> k_mat;
        std::vector<double> f_vec;
        std::vector<std::vector<double>> jac;
        std::vector<double> x_curr;
        std::vector<double> x_new;

        std::vector<double> dt_vec;
        std::vector<double> max_dsw_vec;
        std::vector<double> max_dpr_vec;

        std::vector<double> qo_vec;
        std::vector<double> qw_vec;
        std::vector<double> pinj_vec;

};

double unit_conv = 0.00852702; // units: bar, mD, cP, m, m3/d

SimModel::SimModel():
    t_end(0.),
    max_dsw(0.),
    max_dpr(0.),
    max_dt(0.),
    min_dt(0.),
    conv_tol(0.),
    use_nr(true),
    dt_curr(0.),

    ni(1),
    nj(1),
    dh(1.),
    di(1.),
    dj(1.),
    dz(1.),
    phi(0.01),
    perm(1.),
    p_init(1.),

    bo(1.),
    bw(1.),
    uo(1.),
    uw(1.),
    kr(),

    qwi(1.),
    pwf(1.),
    rw(0.1),
    skin(0.),
    wi(1.),

    t_vec(),
    k_mat(),
    f_vec(),
    jac(),
    x_curr(),
    x_new(),

    dt_vec(),
    max_dsw_vec(),
    max_dpr_vec(),

    qo_vec(),
    qw_vec(),
    pinj_vec()
    {};

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
void SimModel::set_use_nr(bool value){
    use_nr = value;
}
void SimModel::set_ni(int value){
    ni = value;
}
void SimModel::set_nj(int value){
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


int SimModel::get_cell_n(int i, int j){
    return i + 1 + j * ni;
}
double SimModel::get_cell_pr(std::vector<double> x, int c){
    return x[2 * c - 2];
}
double SimModel::get_cell_pr(std::vector<double> x, int i, int j){
    int c = get_cell_n(i,j);
    return get_cell_pr(x, c);
}
double SimModel::get_cell_sw(std::vector<double> x, int c){
    return x[2 * c - 1];
}
double SimModel::get_cell_sw(std::vector<double> x, int i, int j){
    int c = get_cell_n(i,j);
    return get_cell_sw(x, c);
}

double SimModel::get_upstream_sw(std::vector<double> x, int c1, int c2){
    if (c1==c2){
        return get_cell_sw(x,c1);
    } else if (get_cell_pr(x,c1) > get_cell_pr(x,c2)){
        return get_cell_sw(x,c1);
    } else{
        return get_cell_sw(x,c2);
    }
}
double SimModel::get_upstream_sw(std::vector<double> x, int i1, int j1, int i2, int j2){
    int c1 = get_cell_n(i1,j1);
    int c2 = get_cell_n(i2,j2);
    return get_upstream_sw(x, c1, c2);
}

double SimModel::get_kro(std::vector<double> x, int c){
    return get_kro(x, c, c);
}
double SimModel::get_kro(std::vector<double> x, int c1, int c2){
    double swi = get_upstream_sw(x,c1,c2);
    return kr.get_kro(swi);
}
double SimModel::get_kro(std::vector<double> x, int i1, int j1, int i2, int j2){
    int c1 = get_cell_n(i1,j1);
    int c2 = get_cell_n(i2,j2);
    return get_kro(x, c1, c2);
}
double SimModel::get_krw(std::vector<double> x, int c){
    return get_krw(x, c, c);
}
double SimModel::get_krw(std::vector<double> x, int c1, int c2){
    double swi = get_upstream_sw(x,c1,c2);
    return kr.get_krw(swi);
}
double SimModel::get_krw(std::vector<double> x, int i1, int j1, int i2, int j2){
    int c1 = get_cell_n(i1,j1);
    int c2 = get_cell_n(i2,j2);
    return get_krw(x, c1, c2);
}
double SimModel::get_dkro(std::vector<double> x, int c1, int c2){
    double swi = get_upstream_sw(x,c1,c2);
    return kr.get_dkro(swi);
}
double SimModel::get_dkro(std::vector<double> x, int i1, int j1, int i2, int j2){
    int c1 = get_cell_n(i1,j1);
    int c2 = get_cell_n(i2,j2);
    return get_dkro(x, c1, c2);
}
double SimModel::get_dkrw(std::vector<double> x, int c1, int c2){
    double swi = get_upstream_sw(x,c1,c2);
    return kr.get_dkrw(swi);
}
double SimModel::get_dkrw(std::vector<double> x, int i1, int j1, int i2, int j2){
    int c1 = get_cell_n(i1,j1);
    int c2 = get_cell_n(i2,j2);
    return get_dkrw(x, c1, c2);
}

double SimModel::get_tr(int c1, int c2){
    if (abs(c1-c2) == ni){
        return perm * di * dz / dj;
    }
    if (abs(c1-c2) == 1){
        return perm * dj * dz / di;
    }
    return 0.;
}
double SimModel::get_tr(int i1, int j1, int i2, int j2){
    int c1 = get_cell_n(i1,j1);
    int c2 = get_cell_n(i2,j2);
    return get_tr(c1, c2);
}

double SimModel::get_tro(std::vector<double> x, int c1, int c2){
    double tr = get_tr(c1,c2);
    if (tr != 0.){
        tr *= unit_conv * get_kro(x,c1,c2) / (bo * uo);
    }
    return tr;
}
double SimModel::get_tro(std::vector<double> x, int i1, int j1, int i2, int j2){
    int c1 = get_cell_n(i1,j1);
    int c2 = get_cell_n(i2,j2);
    return get_tro(x, c1, c2);
}
double SimModel::get_trw(std::vector<double> x, int c1, int c2){
    double tr = get_tr(c1,c2);
    if (tr != 0.){
        tr *= unit_conv * get_krw(x,c1,c2) / (bw * uw);
    }
    return tr;
}
double SimModel::get_trw(std::vector<double> x, int i1, int j1, int i2, int j2){
    int c1 = get_cell_n(i1,j1);
    int c2 = get_cell_n(i2,j2);
    return get_trw(x, c1, c2);
}

void SimModel::build_K(std::vector<double> x, double dt){
    initial_matrix(k_mat, 2*ni*nj, 2*ni*nj);
    std::vector<int> c_list;
    int nc;
    double tro;
    double trw;
    int c;

    double vp_dt = di * dj * dz * phi / dt;
    double vp_dt_bo = vp_dt / bo;
    double vp_dt_bw = (-1.) * vp_dt / bw;
    for (int j=0; j<nj; j++){
        for (int i=0; i<ni; i++){
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
            nc = static_cast<int>(c_list.size());
            for (int k=0; k<nc; k++){
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
void SimModel::build_K(){
    build_K(x_new, dt_curr);
}
void SimModel::build_F(std::vector<double> x, double dt){
    initial_vector(f_vec, 2*ni*nj);
    int c;

    double vp_dt = di * dj * dz * phi / dt;
    double vp_dt_bo = vp_dt / bo;
    double vp_dt_bw = (-1.) * vp_dt / bw;
    double sw;
    for (int j=0; j<nj; j++){
        for (int i=0; i<ni; i++){
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
void SimModel::build_F(){
    build_F(x_curr, dt_curr);
}
void SimModel::build_Jac(){

}

void SimModel::advance_one_time_step(){

}
void SimModel::save_result(){

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

    double swi = kr.get_swi();
    x_curr.clear();
    for (int j=0; j<nj; j++){
        for (int i=0; i<ni; i++){
            x_curr.push_back(p_init);
            x_curr.push_back(swi);
        }
    }
    x_new = x_curr;
}

// void tests(
//     std::ostream& output,
//     const bool is_2D,
//     const double conv_tol_p,
//     const double conv_tol_sw,
//     const double max_dt,
//     const double min_dt,
//     const size_t max_iterations){

//     std::vector<std::vector<double>> k_matrix;
//     std::vector<double> f_vector;
//     std::vector<double> x_out;

//     std::vector<int> a = {3, 5, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100};
//     for (int i : a) {
//         n = std::to_string(i);

//         output << "Number of cells: " << n << " x " << n2 << std::endl;

//             output << "  Gauss Elimination with Partial Pivot and Scaling" << std::endl;
//             start = high_resolution_clock::now();
//             x_out = SolveGauss(k_matrix, f_vector, true);
//             stop = high_resolution_clock::now();
//             duration = duration_cast<milliseconds>(stop - start);
//             printVector(x_out, output);
//             output << "Elapsed time: " << duration.count() << " milliseconds" << std::endl;
//             check_result(k_matrix, f_vector, x_out, x_true, output);
//             output << std::endl;
//         }

//     }
// }

int main(){
    #ifdef _WIN32
        system("cls");
    #elif __unix__
        system("clear");
    #else
        std::cout << "Running on an unknown system" << std::endl;
    #endif

    // test01();
    std::ofstream outputFile("results_NR_1D.txt");
    if (outputFile.is_open()) {

        outputFile.close();
    } else {
        std::cerr << "Unable to open the file for writing.\n";
    }
}