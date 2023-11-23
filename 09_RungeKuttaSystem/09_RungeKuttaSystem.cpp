#include <iostream>
#include <cmath>
#include <limits>
#include <cstdio>
#include <string>
#include <cstdlib>
#include <vector>
#include <functional>
#include <cstring>
#include <fstream>
#include <iomanip>
using namespace std;

const double eps = 1E-12;

double* copy_array(double *original, int number_points){
    if (number_points < 1){
        return nullptr;
    }
    double* out = new double[number_points];
    for(int i=0; i<number_points; i++){
        out[i] = original[i];
    }
    return out;
}

int* copy_array(int *original, int number_points){
    if (number_points < 1){
        return nullptr;
    }
    int* out = new int[number_points];
    for(int i=0; i<number_points; i++){
        out[i] = original[i];
    }
    return out;
}

double max_array(double *vector, int number_points){
    if (number_points < 1){
        return std::numeric_limits<double>::quiet_NaN();
    }
    double out = vector[0];
    for(int i=1; i<number_points; i++){
        if (vector[i] > out || isnan(out)){
            out = vector[i];
        };
    }
    return out;
}

double* nan_array(int number_points){
    double* out = new double[number_points];
    for (int i = 0; i < number_points; i++) {
        out[i] = std::numeric_limits<double>::quiet_NaN();
    }
    return out;
}

std::vector<double> addVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    std::vector<double> result;
    if (vector1.size() == vector2.size()) {
        result.reserve(vector1.size());
        for (size_t i = 0; i < vector1.size(); ++i) {
            result.push_back(vector1[i] + vector2[i]);
        }
    } else {
        std::cerr << "Error: Vectors have different sizes (" << vector1.size() << " and " << vector2.size() << ") and cannot be added.\n";
    }
    return result;
}

std::vector<double> multiplyVector(const std::vector<double>& vector1, const double alpha) {
    std::vector<double> result;
    result.reserve(vector1.size());
    for (size_t i = 0; i < vector1.size(); ++i) {
        result.push_back(vector1[i] * alpha);
    }
    return result;
}

std::vector<double> weightedSum(const std::vector<std::vector<double>>& vectors, const std::vector<double>& weights, double mutiplier) {
    if (vectors.size() == weights.size() && !vectors.empty()) {
        std::vector<double> result(vectors[0].size(), 0.0);

        for (size_t i = 0; i < vectors.size(); ++i) {
            for (size_t j = 0; j < vectors[i].size(); ++j) {
                result[j] += vectors[i][j] * weights[i];
            }
        }

        for (double& value : result) {
            value *= mutiplier;
        }

        return result;
    } else {
        std::cerr << "Error: Vectors and weights have different sizes or are empty.\n";
        return {};
    }
}
std::vector<double> weightedMean(const std::vector<std::vector<double>>& vectors, const std::vector<double>& weights) {
    double sumOfWeights = 0.0;
    for (double weight : weights) {
        sumOfWeights += weight;
    }
    return weightedSum(vectors, weights, 1./sumOfWeights);
}

char* concatenateStrings(const char* str1, const char* str2) {
    size_t len1 = strlen(str1);
    size_t len2 = strlen(str2);
    char* concatenatedStr = new char[len1 + len2 + 1];
    strcpy(concatenatedStr, str1);
    strcat(concatenatedStr, str2);
    return concatenatedStr;
}

class IVPSystem{
    public:
        IVPSystem();
        using FunctionType = std::function<double(double, const std::vector<double>&)>;
        using ExactType = std::function<double(double)>;

        void set_t_initial(double t);
        void set_t_end(double t);
        void set_time_steps(int n);
        void set_relative_error(bool relative);

        void addFunction(const FunctionType& func);
        void addExact(const ExactType& func);
        void add_y_initial(double value);

        void resetFunctions();
        void resetExactFunctions();
        void reset_y_initial();

        std::vector<double> get_f(double t, const std::vector<double>& values) const;
        std::vector<double> get_exact(double t) const;

        std::vector<double> get_t() const;
        std::vector<std::vector<double>> get_y_aprox() const;
        std::vector<std::vector<double>> get_y_exact() const;
        std::vector<std::vector<double>> get_y_error() const;
        std::vector<int> get_f_evaluations() const;
        int get_current_f_evaluations() const;

        void reset_f_evaluations();
        void initialize_problem();

        void calculate_exact_error();
        void print_solution() const;
        void print_solution(const char* filename) const;

        void solve_euler();
        void solve_rungekutta();
        void solve_rungekutta(int n);
    private:
        void reset_results();
        void reset_error_estimate();
        bool check_vector_sizes();

        std::vector<double> calculate_rungekutta_k(double t, const std::vector<std::vector<double>>& k, std::vector<double> k_multipliers) const;
        void solve_rungekuttaN(bool one_step, const std::vector<double> h_multipliers, const std::vector<std::vector<double>> k_multipliers, const std::vector<double> sum_weights);
        void solve_euler(bool one_step);
        void solve_rungekutta2();
        void solve_rungekutta2(bool one_step);
        void solve_rungekutta3();
        void solve_rungekutta3(bool one_step);
        void solve_rungekutta4();
        void solve_rungekutta4(bool one_step);
        void solve_rungekutta5();
        void solve_rungekutta5(bool one_step);

        bool has_f() const;
        bool has_exact() const;
        bool has_solution() const;
        bool has_exact_error() const;
        bool has_y_initial() const;

        void calculate_exact();

        void print_solution(std::ostream& output) const;

        std::vector<FunctionType> functions;
        std::vector<ExactType> exact_functions;
        std::vector<double> y_initial;

        double t_initial;
        double t_end;
        int time_steps;
        double h_current;
        bool relative_error;

        std::vector<double> t_list;
        std::vector<std::vector<double>> y_aprox;
        std::vector<std::vector<double>> y_exact;
        std::vector<std::vector<double>> y_error;
        std::vector<int> f_evaluations;
};

IVPSystem::IVPSystem():
    functions(),
    exact_functions(),
    y_initial(),
    t_initial(0),
    t_end(1),
    time_steps(100),
    h_current(1.),
    relative_error(false),
    t_list(),
    y_aprox(),
    y_exact(),
    y_error(),
    f_evaluations()
    {};

void IVPSystem::reset_results(){
    t_list.clear();
    y_aprox.clear();
    reset_error_estimate();
    reset_f_evaluations();
}
void IVPSystem::reset_f_evaluations(){
    f_evaluations.clear();
}
void IVPSystem::reset_error_estimate(){
    y_exact.clear();
    y_error.clear();
}
void IVPSystem::initialize_problem(){
    reset_results();
    t_list.push_back(t_initial);
    y_aprox.push_back(y_initial);
    f_evaluations.push_back(0);
    h_current = (t_end - t_initial) / time_steps;
}
bool IVPSystem::check_vector_sizes(){
    bool result = true;
    if (functions.size() != y_initial.size()){
        std::cout << "Number of functions (" << functions.size() << ") does not match the number of initial values (" << y_initial.size() << ")." << std::endl;
        result = false;
    }
    if (has_exact()){
        if (functions.size() != exact_functions.size()){
            std::cout << "Number of functions (" << functions.size() << ") does not match the number of exact functions (" << y_initial.size() << ")." << std::endl;
            result = false;
        }
    }
    if (has_solution()){
        if (t_list.size() != y_aprox.size()){
            std::cout << "Number of functions (" << functions.size() << ") does not match the number of exact functions (" << y_initial.size() << ")." << std::endl;
            result = false;
        }
    }
    return result;
}

void IVPSystem::set_t_initial(double t){
    t_initial = t;
    if (abs(t-t_end)<eps){
        t_end = t_initial + 1;
        std::cout << "Cannot set initial time equal to end time. Setting end time to initial time + 1: " << t_end << std::endl;
    }
    reset_results();
}
void IVPSystem::set_t_end(double t){
    if (abs(t-t_initial)<eps){
        std::cout << "Cannot set end time equal to initial time. Keeping previous value: " << t_end << std::endl;
    } else{
        t_end = t;
    }
    reset_results();
}
void IVPSystem::set_time_steps(int n){
    if (n<1){
        std::cout << "Cannot set less than one time-step. Keeping previous value: " << time_steps << std::endl;
    } else{
        time_steps = n;
    }
    reset_results();
}
void IVPSystem::set_relative_error(bool relative){
    if (relative_error != relative){
        relative_error = relative;
        y_error.clear();
    }
}

void IVPSystem::addFunction(const FunctionType& func){
    functions.push_back(func);
}
void IVPSystem::addExact(const ExactType& func){
    exact_functions.push_back(func);
}
void IVPSystem::add_y_initial(double value){
    y_initial.push_back(value);
}
void IVPSystem::resetFunctions(){
    functions.clear();
}
void IVPSystem::resetExactFunctions(){
    exact_functions.clear();
}
void IVPSystem::reset_y_initial(){
    y_initial.clear();
}

bool IVPSystem::has_f() const{
    return !functions.empty();
}
bool IVPSystem::has_exact() const{
    return !exact_functions.empty();
}
bool IVPSystem::has_solution() const{
    return !y_aprox.empty();
}
bool IVPSystem::has_exact_error() const{
    return !y_error.empty();
}
bool IVPSystem::has_y_initial() const{
    return !y_initial.empty();
}

std::vector<double> IVPSystem::get_f(double t, const std::vector<double>& values) const{
    std::vector<double> results;
    for (const auto& func: functions) {
        double result = func(t, values);
        results.push_back(result);
    }
    return results;
}
std::vector<double> IVPSystem::get_exact(double t) const{
    std::vector<double> results;
    for (const auto& func: exact_functions) {
        double result = func(t);
        results.push_back(result);
    }
    return results;
}

std::vector<double> IVPSystem::get_t() const{
    return t_list;
}
std::vector<std::vector<double>> IVPSystem::get_y_aprox() const{
    return y_aprox;
}
std::vector<std::vector<double>> IVPSystem::get_y_exact() const{
    return y_exact;
}
std::vector<std::vector<double>> IVPSystem::get_y_error() const{
    return y_error;
}
std::vector<int> IVPSystem::get_f_evaluations() const{
    return f_evaluations;
}

int IVPSystem::get_current_f_evaluations() const{
    if (!f_evaluations.empty()) {
        return f_evaluations.back();
    } else {
        return 0;
    }
}

void IVPSystem::calculate_exact(){
    if (!has_exact()){
        std::cout << "Exact solution not defined. Cannot continue." << std::endl;
        return;
    }
    y_exact.clear();
    for (const auto& t: t_list){
        y_exact.push_back(get_exact(t));
    }
}

void IVPSystem::calculate_exact_error(){
    if (!has_solution()){
        std::cout << "No solutions found. Cannot continue." << std::endl;
        return;
    }
    if (!has_exact()){
        std::cout << "Exact solution not defined. Cannot continue." << std::endl;
        return;
    }
    if (!check_vector_sizes()){
        return;
    }

    calculate_exact();
    y_error.clear();
    std::vector<double> error_list;
    double error_single;
    for (size_t i = 0; i < y_aprox.size(); ++i){
        error_list.clear();
        for (size_t j = 0; j < y_aprox[i].size(); ++j){
            error_single = abs(y_aprox[i][j] - y_exact[i][j]);
            if (relative_error){
                if (abs(y_exact[i][j]) > eps){
                    error_single = error_single / abs(y_exact[i][j]);
                } else if(error_single > eps){
                    error_single = std::numeric_limits<double>::quiet_NaN();
                } else{
                    error_single = 0.;
                }
            }
            error_list.push_back(error_single);
        }
        y_error.push_back(error_list);
    }
}

void IVPSystem::print_solution() const{
    print_solution(std::cout);
}

void IVPSystem::print_solution(const char* filename) const{
    std::ofstream outputFile(filename);
    if (outputFile.is_open()) {
        print_solution(outputFile);
        outputFile.close();
    } else {
        std::cerr << "Unable to open the file for writing.\n";
    }
}

void IVPSystem::print_solution(std::ostream& output) const{
    if (!has_solution()){
        std::cout << "There is no solution to print." << std::endl;
        return;
    }

    output << std::setw(3) << "i";
    output << "\t" << std::setw(16) << "t";
    output << "\t" << std::setw(16) << "evaluations";
    for (size_t i = 1; i <= functions.size(); ++i) {
        output << "\t" << std::setw(16) << "y_approx" << i;
    }
    if (has_exact_error()){
        for (size_t i = 1; i <= exact_functions.size(); ++i) {
            output << "\t" << std::setw(16) << "y_exact" << i;
        }
        for (size_t i = 1; i <= exact_functions.size(); ++i) {
            output << "\t" << std::setw(16) << "y_error" << i;
        }
    }
    output << "\n";

    for (size_t i = 0; i < t_list.size(); ++i) {
        output << std::setw(3) << i;
        output << "\t" << std::setw(16) << t_list[i];
        output << "\t" << std::setw(16) << f_evaluations[i];
        for (double value: y_aprox[i]) {
            output << "\t" << std::setw(16) << value;
        }
        if (has_exact_error()){
            for (double value: y_exact[i]) {
                output << "\t" << std::setw(16) << value;
            }
            for (double value: y_error[i]) {
                output << "\t" << std::setw(16) << value;
            }
        }
        output << "\n";
    }
}

void IVPSystem::solve_rungekutta(){
    solve_rungekutta(4);
}
void IVPSystem::solve_rungekutta(int n){
    initialize_problem();
    switch (n)
    {
    case 1:
        solve_euler();
        break;
    case 2:
        solve_rungekutta2();
        break;
    case 3:
        solve_rungekutta3();
        break;
    case 4:
        solve_rungekutta4();
        break;
    case 5:
        solve_rungekutta5();
        break;
    default:
        std::cout << "No Runge-Kutta method with order " << n << " was implemented.\n";
        break;
    }
}

void IVPSystem::solve_euler(){
    initialize_problem();
    solve_euler(false);
}
void IVPSystem::solve_euler(bool one_step){
    solve_rungekuttaN(one_step, {1.}, {{}}, {1.});
}
void IVPSystem::solve_rungekutta2(){
    initialize_problem();
    solve_rungekutta2(false);
}
void IVPSystem::solve_rungekutta2(bool one_step){
    solve_rungekuttaN(one_step, {0, 0.5}, {{},{0.5}}, {0., 1.});
}
void IVPSystem::solve_rungekutta3(){
    initialize_problem();
    solve_rungekutta3(false);
}
void IVPSystem::solve_rungekutta3(bool one_step){
    solve_rungekuttaN(one_step, {0, 1./3., 2./3.}, {{},{1./3.},{0.,2./3.}}, {0.25, 0., 0.75});
}
void IVPSystem::solve_rungekutta4(){
    initialize_problem();
    solve_rungekutta4(false);
}
void IVPSystem::solve_rungekutta4(bool one_step){
    solve_rungekuttaN(one_step, {0, 0.5, 0.5, 1}, {{},{0.5},{0.,0.5},{0.,0.,1.}}, {1., 2., 2., 1.});
}
void IVPSystem::solve_rungekutta5(){
    initialize_problem();
    solve_rungekutta5(false);
}
void IVPSystem::solve_rungekutta5(bool one_step){
    std::vector<std::vector<double>> k_multiplier;
    k_multiplier.clear();
    k_multiplier.push_back({});
    k_multiplier.push_back({0.25});
    k_multiplier.push_back({3./32., 9./32.});
    k_multiplier.push_back({1932./2197., -7200./2197., 7296./2197.});
    k_multiplier.push_back({439./216., -8., 3680./513., -845./4104.});
    k_multiplier.push_back({-8./27., 2., -3544./2565., 1859./4104., -11./40.});
    solve_rungekuttaN(one_step, {0, 0.25, 3./8., 12./13., 1., 0.5}, k_multiplier, {16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.});
}

std::vector<double> IVPSystem::calculate_rungekutta_k(double t, const std::vector<std::vector<double>>& k, std::vector<double> k_multipliers) const{
    std::vector<double> results, values, k_in;
    values = y_aprox.back();
    for (size_t i = 0; i < k_multipliers.size(); ++i){
        k_in = multiplyVector(k[i], k_multipliers[i]);
        values = addVectors(values, k_in);
    }
    for (const auto& func: functions) {
        double result = h_current * func(t, values);
        results.push_back(result);
    }
    return results;
}

void IVPSystem::solve_rungekuttaN(bool one_step, const std::vector<double> h_multipliers, const std::vector<std::vector<double>> k_multipliers, const std::vector<double> sum_weights){
    if (!has_f()){
        std::cout << "Functions f(t,y) not defined. Cannot continue." << std::endl;
        return;
    }

    int current_step = static_cast<int>(t_list.size()) - 1;
    if (current_step < 0){
        initialize_problem();
        current_step++;
    }
    int steps = time_steps;
    if (one_step){
        steps = current_step+1;
    }
    std::vector<std::vector<double>> k;
    std::vector<double> w;
    double t, t_rk;
    for (int i=current_step; i<steps; i++){
        t = t_initial + i*h_current;
        k.clear();
        for (size_t j = 0; j < h_multipliers.size(); ++j){
            t_rk = t + h_current * h_multipliers[j];
            k.push_back(calculate_rungekutta_k(t_rk, k, k_multipliers[j]));
        }
        w = weightedMean(k, sum_weights);
        w = addVectors(y_aprox.back(), w);

        f_evaluations.push_back(f_evaluations.back() + static_cast<int>(h_multipliers.size()));
        t_list.push_back(t+h_current);
        y_aprox.push_back(w);
    }
    reset_error_estimate();
}


// ############# Tests #############

double f_cumulative(double t, std::vector<double> y){
    return y[0] + 0*t;
}

double f_test_1(double t, std::vector<double> y){
    return y[0] - t*t + 1;
}
double exact_test_1(double t){
    return pow(t+1, 2) - 0.5 * exp(t); // y(0)=0.5
}
double exact_test_1_cum(double t){
    return t*(1. + t*(1. + t/3.) ) + 0.5 *(1. - exp(t)); //y(0)=0
}

double f_test_2(double t, std::vector<double> y){
    return -2.*y[0] + 3.*exp(t);
}
double exact_test_2(double t){
    return 2. * exp(-2.*t) + exp(t);  // y(0)=3
}
double exact_test_2_cum(double t){
    return - exp(-2.*t) + exp(t);  // y(0)=0
}

double f_test_3(double t, std::vector<double> y){
    return 4*cos(t) - 8*sin(t) + 2*y[0];
}
double exact_test_3(double t){
    return 4*sin(t) + 3*exp(2*t);  // y(0)=3
}
double exact_test_3_cum(double t){
    return 4 - 4*cos(t) + 3*exp(t)*sinh(t); //y(0)=0
}

double a = 3.;
double f_test_4(double t, std::vector<double> y){
    return -a*y[0] + 0*t;
}
double exact_test_4(double t){
    return 3.*exp(-a*t);  // y(0)=3
}
double exact_test_4_cum(double t){
    return 3/a * (1- exp(-a*t));  // y(0) = 0
}

double f_test_5A(double t, std::vector<double> y){
    return -4.*y[0] + 3.*y[1] + 6. + 0*t;
}
double f_test_5B(double t, std::vector<double> y){
    return -2.4*y[0] + 1.6*y[1] + 3.6 + 0*t;
}
double exact_test_5A(double t){
    return -3.375*exp(-2.*t) + 1.875*exp(-0.4*t) + 1.5;  // y(0)=0
}
double exact_test_5B(double t){
    return -2.25*exp(-2.*t) + 2.25*exp(-0.4*t);  // y(0)=0
}

double f_test_6A(double t, std::vector<double> y){
    return y[1] + 0.*t;
}
double f_test_6B(double t, std::vector<double> y){
    return exp(2.*t) * sin(t) - 2.*y[0] + 2.*y[1];
}
double exact_test_6A(double t){
    return 0.2*exp(2.*t) * (sin(t) - 2.*cos(t));  // y(0)=-0.4
}
double exact_test_6B(double t){
    return 0.2*exp(2.*t) * (4.*sin(t) - 3.*cos(t));  // y(0)=-0.6
}


void test_rungekutta(IVPSystem ivp, string problemname, const char* filename, FILE* evaluationsFile){
    const char* name = problemname.c_str();
    ivp.set_relative_error(true);

    const int n = 10;

    printf(" Problem %s: Runge-Kutta 4\n", name);
    ivp.set_time_steps(6*n);
    ivp.reset_f_evaluations();
    ivp.solve_rungekutta(4);
    printf("    'f' evaluations: %d\n", ivp.get_current_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "RungeKutta4", ivp.get_current_f_evaluations());
    ivp.calculate_exact_error();
    ivp.print_solution(concatenateStrings(filename,"_rungekutta4.txt"));

    printf(" Problem %s: Runge-Kutta 5\n", name);
    ivp.set_time_steps(4*n);
    ivp.reset_f_evaluations();
    ivp.solve_rungekutta(5);
    printf("    'f' evaluations: %d\n", ivp.get_current_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "RungeKutta5", ivp.get_current_f_evaluations());
    ivp.calculate_exact_error();
    ivp.print_solution(concatenateStrings(filename,"_rungekutta5.txt"));

    printf(" Problem %s: Runge-Kutta 5 'Best'\n", name);
    ivp.set_time_steps(100*n);
    ivp.reset_f_evaluations();
    ivp.solve_rungekutta(5);
    printf("    'f' evaluations: %d\n", ivp.get_current_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "RungeKutta5Best", ivp.get_current_f_evaluations());
    ivp.calculate_exact_error();
    ivp.print_solution(concatenateStrings(filename,"_rungekutta5Best.txt"));
}

void tests_rungekutta(){

    FILE* outFile = fopen("evaluations.txt", "w");
    if (outFile == nullptr) {
        printf("Coud not create file: %s\n","evaluations.txt");
    }
    fprintf(outFile, "%s\t%s\t%s\n", "Problem", "Method", "Evaluations");

    printf("### Problem #1 ###\n");
    IVPSystem test1;
    test1.addFunction(f_test_1);
    test1.add_y_initial(0.5);
    test1.addFunction(f_cumulative);
    test1.add_y_initial(0.0);
    test1.set_t_initial(0.);
    test1.set_t_end(2.);
    test1.addExact(exact_test_1);
    test1.addExact(exact_test_1_cum);
    test_rungekutta(test1, "#1", "test1", outFile);

    printf("### Problem #2 ###\n");
    IVPSystem test2;
    test2.addFunction(f_test_2);
    test2.add_y_initial(3.);
    test2.addFunction(f_cumulative);
    test2.add_y_initial(0.);
    test2.set_t_initial(0.);
    test2.set_t_end(2.);
    test2.addExact(exact_test_2);
    test2.addExact(exact_test_2_cum);
    test_rungekutta(test2, "#2", "test2", outFile);

    printf("### Problem #3 ###\n");
    IVPSystem test3;
    test3.addFunction(f_test_3);
    test3.add_y_initial(3.);
    test3.addFunction(f_cumulative);
    test3.add_y_initial(0.);
    test3.set_t_initial(0.);
    test3.set_t_end(2.);
    test3.addExact(exact_test_3);
    test3.addExact(exact_test_3_cum);
    test_rungekutta(test3, "#3", "test3", outFile);

    printf("### Problem #4 ###\n");
    IVPSystem test4;
    test4.addFunction(f_test_4);
    test4.add_y_initial(3.);
    test4.addFunction(f_cumulative);
    test4.add_y_initial(0.);
    test4.set_t_initial(0.);
    test4.set_t_end(2.);
    test4.addExact(exact_test_4);
    test4.addExact(exact_test_4_cum);
    test_rungekutta(test4, "#4", "test4", outFile);

    printf("### Problem #5 ###\n");
    IVPSystem test5;
    test5.addFunction(f_test_5A);
    test5.add_y_initial(0.);
    test5.addFunction(f_test_5B);
    test5.add_y_initial(0.);
    test5.set_t_initial(0.);
    test5.set_t_end(2.);
    test5.addExact(exact_test_5A);
    test5.addExact(exact_test_5B);
    test_rungekutta(test5, "#5", "test5", outFile);

    printf("### Problem #6 ###\n");
    IVPSystem test6;
    test6.addFunction(f_test_6A);
    test6.add_y_initial(-0.4);
    test6.addFunction(f_test_6B);
    test6.add_y_initial(-0.6);
    test6.set_t_initial(0.);
    test6.set_t_end(2.);
    test6.addExact(exact_test_6A);
    test6.addExact(exact_test_6B);
    test_rungekutta(test6, "#6", "test6", outFile);

    fclose(outFile);
}

class Fetkovich{
    public:
        Fetkovich();
        void set_aquifer_productivity_index(double value); // m3/d / bar
        void set_aquifer_initial_pressure(double value);   // bar
        void set_aquifer_initial_pore_volume(double value);  // m3
        void set_aquifer_total_compressibility(double value);  // 1/bar

        void set_reservoir_initial_pressure(double value);   // bar
        void set_reservoir_pressure_function(double (*f)(double, double)); // f(We, t)
        void set_exact_water_flow_function(double (*f)(double)); // f(t)
        void set_exact_water_cumulative_function(double (*f)(double)); // f(t)

        double get_aquifer_flow_rate(double t, double pr) const;  // m3/d
        double get_aquifer_delta_cumulative_flow(double dt, double paq_avg, double pr_avg) const;  // m3

        void solve_aquifer_flow(double t_end, int steps);
        double* get_result_times() const; // d
        double* get_result_water_flow() const; // m3/d
        double* get_result_cumulative_flow() const; // m3
        double* get_result_aquifer_pressure() const; // bar
        double* get_result_reservoir_pressure() const; // bar

        double* get_result_water_flow_error(bool relative) const; // m3/d or adm.
        double* get_result_cumulative_flow_error(bool relative) const; // m3 or adm.

        double get_exact(double t) const;
        double get_exact_cumulative(double t) const;

        void print_solution();
        void print_solution(string filename);

        int get_f_evaluations() const;
        void reset_f_evaluations();

    private:
        double We_max() const;

        double get_aquifer_pressure(double We) const;   // bar

        bool has_f_pr() const;
        double get_f_pr(double We, double t) const;

        bool has_exact() const;
        bool has_exact_cumulative() const;

        bool has_solution() const;

        double J;
        double pi;
        double Wi;
        double ct;
        double pr_0;

        double (*f_pr_pointer)(double, double);

        int time_steps;
        double *t_pointer;
        double *Qw_pointer;
        double *We_pointer;
        double *p_aquifer_pointer;
        double *p_reservoir_pointer;

        double (*f_exact_pointer)(double);
        double (*f_exact_cum_pointer)(double);

        int f_evaluations;
};

Fetkovich::Fetkovich():
    J(0.), pi(0.), Wi(0.), ct(0.), pr_0(0.),
    f_pr_pointer(nullptr),
    time_steps(0),
    t_pointer(nullptr),
    Qw_pointer(nullptr),
    We_pointer(nullptr),
    p_aquifer_pointer(nullptr),
    p_reservoir_pointer(nullptr),
    f_exact_pointer(nullptr),
    f_exact_cum_pointer(nullptr),
    f_evaluations(0)
    {};

void Fetkovich::set_aquifer_productivity_index(double value){
    J = value;
}
void Fetkovich::set_aquifer_initial_pressure(double value){
    pi = value;
}
void Fetkovich::set_aquifer_initial_pore_volume(double value){
    Wi = value;
}
void Fetkovich::set_aquifer_total_compressibility(double value){
    ct = value;
}
void Fetkovich::set_reservoir_initial_pressure(double value){
    pr_0 = value;
}
double Fetkovich::We_max() const{
    return Wi * ct * pi;
}

double Fetkovich::get_aquifer_flow_rate(double t, double pr) const{
    return exp(-J * t * pi / We_max()) * J * (pi - pr);
}
double Fetkovich::get_aquifer_delta_cumulative_flow(double dt, double paq_avg, double pr_avg) const{
    return (1 - exp(-J * dt * pi/ We_max())) * ct * Wi * (paq_avg - pr_avg);
}
double Fetkovich::get_aquifer_pressure(double We) const{
    return pi *( 1 - We / We_max());
}

void Fetkovich::set_reservoir_pressure_function(double (*f)(double, double)){
    f_pr_pointer = f;
}
bool Fetkovich::has_f_pr() const{
    return f_pr_pointer != nullptr;
}
double Fetkovich::get_f_pr(double t, double We) const{
    if (!has_f_pr()){
        printf("Reservoir pressure function not defined.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return f_pr_pointer(t,We);
}

void Fetkovich::set_exact_water_flow_function(double (*f)(double)){
    f_exact_pointer = f;
}
bool Fetkovich::has_exact() const{
    return f_exact_pointer != nullptr;
}
double Fetkovich::get_exact(double t) const{
    if (!has_exact()){
        printf("Exact function not defined.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return f_exact_pointer(t);
}

void Fetkovich::set_exact_water_cumulative_function(double (*f)(double)){
    f_exact_cum_pointer = f;
}
bool Fetkovich::has_exact_cumulative() const{
    return f_exact_pointer != nullptr;
}
double Fetkovich::get_exact_cumulative(double t) const{
    if (!has_exact_cumulative()){
        printf("Exact cumulative function not defined.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return f_exact_cum_pointer(t);
}

int Fetkovich::get_f_evaluations() const{
    return f_evaluations;
}
void Fetkovich::reset_f_evaluations(){
    f_evaluations = 0;
}

bool Fetkovich::has_solution() const{
    return t_pointer != nullptr;
}
void Fetkovich::solve_aquifer_flow(double t_end, int steps){
    time_steps = steps+1;
    double* t = new double[steps+1];
    double* Qw = new double[steps+1];
    double* We = new double[steps+1];
    double* p_aq = new double[steps+1];
    double* p_res = new double[steps+1];

    t[0] = 0;
    Qw[0] = J * (pi - pr_0);
    We[0] = 0;
    p_aq[0] = pi;
    p_res[0] = pr_0;

    double dWe = 0;
    double pr_avg = pr_0;
    double pr_avg_old;
    double count;
    double dt = t_end / steps;

    for (int i = 1; i < steps+1; i++) {
        t[i] = i * dt;

        pr_avg_old = 0;
        count = 0;
        while (abs(pr_avg_old - pr_avg) > 1E-1 && count<=20){
            pr_avg_old = pr_avg;
            dWe = get_aquifer_delta_cumulative_flow(dt, p_aq[i-1], pr_avg);
            f_evaluations++;
            We[i] = We[i-1] + dWe;
            p_res[i] = get_f_pr(t[i], We[i]);
            pr_avg = (p_res[i] + p_res[i-1])/2;
            count++;
        }
        p_aq[i] = get_aquifer_pressure(We[i]);
        Qw[i] = J * (p_aq[i] - p_res[i]);
    }

    t_pointer = t;
    Qw_pointer = Qw;
    We_pointer = We;
    p_aquifer_pointer = p_aq;
    p_reservoir_pointer = p_res;
}
double* Fetkovich::get_result_times() const{
    return t_pointer;
}
double* Fetkovich::get_result_water_flow() const{
    return Qw_pointer;
}
double* Fetkovich::get_result_cumulative_flow() const{
    return We_pointer;
}
double* Fetkovich::get_result_aquifer_pressure() const{
    return p_aquifer_pointer;
}
double* Fetkovich::get_result_reservoir_pressure() const{
    return p_reservoir_pointer;
}

double* Fetkovich::get_result_water_flow_error(bool relative) const{
    if (!has_solution()){
        printf("There is no solution. Cannot continue.\n");
        return nan_array(time_steps);
    }
    if (!has_exact()){
        printf("There is no exact solution. Cannot continue.\n");
        return nan_array(time_steps);
    }

    double exact, error;
    double* error_list = new double[time_steps];
    for (int i=0; i<time_steps; i++){
        exact = get_exact(t_pointer[i]);
        error = abs(exact - Qw_pointer[i]);
        if (relative){
            if (abs(exact) > eps){
                error = error/abs(exact);
            } else{
                error = std::numeric_limits<double>::quiet_NaN();
            }
        }
        error_list[i] = error;
    }
    return error_list;
}

double* Fetkovich::get_result_cumulative_flow_error(bool relative) const{
    if (!has_solution()){
        printf("There is no solution. Cannot continue.\n");
        return nan_array(time_steps);
    }
    if (!has_exact_cumulative()){
        printf("There is no exact solution. Cannot continue.\n");
        return nan_array(time_steps);
    }

    double exact, error;
    double* error_list = new double[time_steps];
    for (int i=0; i<time_steps; i++){
        exact = get_exact_cumulative(t_pointer[i]);
        error = abs(exact - We_pointer[i]);
        if (relative){
            if (abs(exact) > eps){
                error = error/abs(exact);
            } else{
                error = std::numeric_limits<double>::quiet_NaN();
            }
        }
        error_list[i] = error;
    }
    return error_list;
}

void Fetkovich::print_solution(){
    if (!has_solution()){
        printf("There is no solution to print.\n");
        return;
    }

    printf("%5s\t%16s","i","Time");
    printf("\t%16s","Res.Pres.");
    printf("\t%16s","Aq.Pres.");
    printf("\t%16s","Wat.Flow");
    printf("\t%16s","Cumulative");
    if (has_exact()){
        printf("\t%16s","ExactQw");
        printf("\t%16s","Error");
        printf("\t%16s","ErrorRel");
    }
    if (has_exact_cumulative()){
        printf("\t%16s","ExactWe");
        printf("\t%16s","ErrorWe");
        printf("\t%16s","ErrorRelWe");
    }
    printf("\n");

    double exact;
    double error;
    for (int i=0; i<time_steps; i++){
        printf("%5i\t%16.10g", i, t_pointer[i]);
        printf("\t%16.10g", p_reservoir_pointer[i]);
        printf("\t%16.10g", p_aquifer_pointer[i]);
        printf("\t%16.10g", Qw_pointer[i]);
        printf("\t%16.10g", We_pointer[i]);
        if (has_exact()){
            exact = get_exact(t_pointer[i]);
            printf("\t%16.10g", exact);
            error = abs(exact - Qw_pointer[i]);
            printf("\t%16.10g", error);
            if (abs(exact) > eps){
                printf("\t%16.10g", error/abs(exact));
            } else{
                printf("\t%16.10g", std::numeric_limits<double>::quiet_NaN());
            }
        }
        if (has_exact_cumulative()){
            exact = get_exact_cumulative(t_pointer[i]);
            printf("\t%16.10g", exact);
            error = abs(exact - We_pointer[i]);
            printf("\t%16.10g", error);
            if (abs(exact) > eps){
                printf("\t%16.10g", error/abs(exact));
            } else{
                printf("\t%16.10g", std::numeric_limits<double>::quiet_NaN());
            }
        }
        printf("\n");
    }
}

void Fetkovich::print_solution(string filename){
    if (!has_solution()){
        printf("There is no solution to print.\n");
        return;
    }

    FILE* outFile = fopen(filename.c_str(), "w");
    if (outFile == nullptr) {
        printf("Coud not create file: %s\n",filename.c_str());
        return;
    }

    fprintf(outFile,"%5s\t%16s","i","Time");
    fprintf(outFile,"\t%16s","Res.Pres.");
    fprintf(outFile,"\t%16s","Aq.Pres.");
    fprintf(outFile,"\t%16s","Wat.Flow");
    fprintf(outFile,"\t%16s","Cumulative");
    if (has_exact()){
        fprintf(outFile,"\t%16s","ExactQw");
        fprintf(outFile,"\t%16s","Error");
        fprintf(outFile,"\t%16s","ErrorRel");
    }
    if (has_exact_cumulative()){
        fprintf(outFile,"\t%16s","ExactWe");
        fprintf(outFile,"\t%16s","ErrorWe");
        fprintf(outFile,"\t%16s","ErrorRelWe");
    }
    fprintf(outFile,"\n");

    double exact;
    double error;
    for (int i=0; i<time_steps; i++){
        fprintf(outFile,"%5i\t%16.10g", i, t_pointer[i]);
        fprintf(outFile,"\t%16.10g", p_reservoir_pointer[i]);
        fprintf(outFile,"\t%16.10g", p_aquifer_pointer[i]);
        fprintf(outFile,"\t%16.10g", Qw_pointer[i]);
        fprintf(outFile,"\t%16.10g", We_pointer[i]);
        if (has_exact()){
            exact = get_exact(t_pointer[i]);
            fprintf(outFile,"\t%16.10g", exact);
            error = abs(exact - Qw_pointer[i]);
            fprintf(outFile,"\t%16.10g", error);
            if (abs(exact) > eps){
                fprintf(outFile,"\t%16.10g", error/abs(exact));
            } else{
                fprintf(outFile,"\t%16.10g", std::numeric_limits<double>::quiet_NaN());
            }
        }
        if (has_exact_cumulative()){
            exact = get_exact_cumulative(t_pointer[i]);
            fprintf(outFile,"\t%16.10g", exact);
            error = abs(exact - We_pointer[i]);
            fprintf(outFile,"\t%16.10g", error);
            if (abs(exact) > eps){
                fprintf(outFile,"\t%16.10g", error/abs(exact));
            } else{
                fprintf(outFile,"\t%16.10g", std::numeric_limits<double>::quiet_NaN());
            }
        }
        fprintf(outFile,"\n");
    }
    fclose(outFile);
}

double Newman_Consolidated_Sandstone(double por){
    if (por >= 1.){
        por = por / 100.;
    }
    return 97.32E-6 / (1 + 55.8721 * pow(por, 1.428586) ) * 14.50377; // 1/bar
}
double Newman_Limestone(double por){
    if (por >= 1.){
        por = por / 100.;
    }
    return 0.8535 / (1 + 2.367E6 * pow(por, 0.93023) ) * 14.50377; // 1/bar
}

double Standing_p_bubble(double api, double dg, double gor, double t){
    double x = 0.0125 * api - 0.00091 * (1.8 * t + 32);
    double y = gor / dg / 0.1373;
    return pow(y, 1/1.205) * pow(10, -x);
}
double Standing_co_bubble(double api, double dg, double gor, double t){
    double x = -1433. + 5 * 5.615 * gor;
    double y = 17.2*(1.8 * t + 491.67);
    double z = -1180. * dg + 12.61 * api;
    double w = 1E5 * Standing_p_bubble(api, dg, gor, t);
    return (x + y + z) / w;
}
double Standing_bo_bubble(double api, double dg, double gor, double t){
    double d_o = 141.5 / (131.5 + api);
    double x = 5.615 * gor * pow(dg / d_o, 0.5);
    double y = 1.25 * (1.8 * t + 32);
    return 0.9759 +  12E-5 * pow(x + y, 1.2);
}


double qo_problem_1(double t){
    if (t < 30.){
        return 500.;
    } else{
        return 0.;
    }
}
double np_problem_1(double t){
    if (t < 30.){
        return 500.*t;
    } else{
        return 500.*30.;
    }
}

double reservoir_problem_1(int var_number, double t, std::vector<double> parameters){
    // Common parameters: oil reservoir
    double api = 25.;
    double dg = 0.6;
    double rgo = 60.;
    double temp = 65.;
    double Bob = Standing_bo_bubble(api, dg, rgo, temp);
    double Cob = Standing_co_bubble(api, dg, rgo, temp);
    double Pb  = Standing_p_bubble(api, dg, rgo, temp);
    double Cpor = Newman_Consolidated_Sandstone(0.2);
    double Bw = 1.;
    double Voil = 0.8 * 1*1E6;
    double Vwat = 0.2 * 1*1E6;
    double p_ini = 230.;

    double VoilIP_0 = Voil * Bob * (1 + Cob*(Pb - p_ini));
    double VwatIP_0 = Vwat * Bw;
    double Vpor_0 = VoilIP_0 + VwatIP_0;

    // Common parameters: aquifer
    double J = 20.;
    double ct = Newman_Consolidated_Sandstone(0.03);
    double Wi = 1.*1E6;
    double paq_ini = p_ini;
    double We_max = Wi * ct * paq_ini;

    double Qw = parameters[0]; // Qw (water flow from aquifer to reservoir)
    double We = parameters[1]; // We (cumulative Qw)
    // double Pr = parameters[2]; // Pr (pressure at the interface)
    double Qo = qo_problem_1(t); // Qo (oil flow from reservoir)
    double Np = np_problem_1(t); // Np (cumulative Qo)

    double beta = (Voil - Np) * Bob * (1. + Cob * Pb) + (Vwat + We) * Bw - Vpor_0 * (1. - Cpor * p_ini);
    double alpha = (Voil - Np) * Bob * Cob + Vpor_0 * Cpor;
    double beta_p = -Qo * Bob * (1. + Cob * Pb) + Qw * Bw;
    double alpha_p = -Qo * Bob * Cob;

    double d_pres_dt = 1./alpha*(beta_p - beta/alpha * alpha_p);

    // Initial conditions
    if (t<0){
        if (var_number == 0){
            return J * (paq_ini - p_ini);
        } if (var_number == 1){
            return 0.;
        }if (var_number == 2){
            return p_ini;
        }
    }

    if (var_number == -1){
        // interface pressure for Fetkovich method
        return beta / alpha;
    } else if (var_number == 0){
        return -J * paq_ini / We_max * Qw - J * d_pres_dt;
    } if (var_number == 1){
        return Qw;
    }if (var_number == 2){
        return d_pres_dt;
    } else{
        return std::numeric_limits<double>::quiet_NaN();
    }
}
double f_pres_problem_1(double t, double We){
    return reservoir_problem_1(-1, t, {0., We, 0., 0.});
}
double f_pres_problem_1_Qw(double t, std::vector<double> parameters){
    return reservoir_problem_1(0, t, parameters);
}
double f_pres_problem_1_We(double t, std::vector<double> parameters){
    return reservoir_problem_1(1, t, parameters);
}
double f_pres_problem_1_Pr(double t, std::vector<double> parameters){
    return reservoir_problem_1(2, t, parameters);
}


void Fetkovich_tests(){
    FILE* outFile = fopen("evaluationsAquifer.txt", "w");
    if (outFile == nullptr) {
        printf("Coud not create file: %s\n","evaluationsAquifer.txt");
    }
    fprintf(outFile, "%s\t%s\t%s\n", "Problem", "Method", "Evaluations");

    printf("Instant Pressure Equilibrium Reservoir with Fetkovich\n");
    Fetkovich aqFet;
    aqFet.set_aquifer_initial_pore_volume(1*1E6);
    aqFet.set_aquifer_initial_pressure(230.);
    aqFet.set_aquifer_productivity_index(20.);
    aqFet.set_aquifer_total_compressibility(Newman_Consolidated_Sandstone(0.03));
    aqFet.set_reservoir_initial_pressure(230.);
    aqFet.set_reservoir_pressure_function(f_pres_problem_1);

    aqFet.solve_aquifer_flow(200., 28*14/2);
    printf("    'f' evaluations: %d\n", aqFet.get_f_evaluations());
    aqFet.print_solution("aq1_fetkovich.txt");

    fprintf(outFile, "%s\t%s\t%d\n", "Aquifer#1", "Fetkovich", aqFet.get_f_evaluations());

    printf("Instant Pressure Equilibrium Reservoir as an IVP\n");
    IVPSystem aqIVP1;
    aqIVP1.addFunction(f_pres_problem_1_Qw);
    aqIVP1.add_y_initial(f_pres_problem_1_Qw(-1.,{0.,0.,0.}));
    aqIVP1.addFunction(f_pres_problem_1_We);
    aqIVP1.add_y_initial(f_pres_problem_1_We(-1.,{0.,0.,0.}));
    aqIVP1.addFunction(f_pres_problem_1_Pr);
    aqIVP1.add_y_initial(f_pres_problem_1_Pr(-1.,{0.,0.,0.}));

    aqIVP1.set_t_initial(0.);
    aqIVP1.set_t_end(200.);

    test_rungekutta(aqIVP1, "Aquifer#1", "aq1", outFile);
    fclose(outFile);

    return;

    // int n_tests = 11;
    // int* steps = new int[n_tests]{5, 10, 20, 50, 80, 100, 150, 200, 250, 500, 1000};

    // double* error_list = new double;
    // double error_end, error_max;
    // int evaluations;

    // printf("We Error Sensibility with Fetkovich\n");
    // outFile = fopen("aq1_fetkovich_sens.txt", "w");
    // // printf("%10s\t%16s\t%16s\t%16s\n","Steps", "Evaluations", "ErrorEnd", "ErrorMax");
    // fprintf(outFile,"%10s\t%16s\t%16s\t%16s\n","Steps", "Evaluations", "ErrorEnd", "ErrorMax");
    // for (int i=0; i<10; i++){
    //     aqFet.reset_f_evaluations();
    //     aqFet.solve_aquifer_flow(200., steps[i]);
    //     evaluations = aqFet.get_f_evaluations();
    //     error_list = aqFet.get_result_cumulative_flow_error(true);
    //     error_end = error_list[steps[i]];
    //     error_max = max_array(error_list, steps[i]+1);
    //     // printf("%10d\t%16d\t%16.10g\t%16.10g\n",steps[i], evaluations, error_end, error_max);
    //     fprintf(outFile,"%10d\t%16d\t%16.10g\t%16.10g\n",steps[i], evaluations, error_end, error_max);
    // }
    // fclose(outFile);

    // string stringArray[2]={"rungekutta4","rungekutta5"};

    // aqIVP1.set_relative_error(true);
    // for (int j=0; j<2; j++){
    //     printf("We Error Sensibility with %s\n", stringArray[j].c_str());
    //     outFile = fopen(("aq1_" + stringArray[j] + "_sens.txt").c_str(), "w");
    //     // printf("%10s\t%16s\t%16s\t%16s\n","Steps", "Evaluations", "ErrorEnd", "ErrorMax");
    //     fprintf(outFile,"%10s\t%16s\t%16s\t%16s\n","Steps", "Evaluations", "ErrorEnd", "ErrorMax");
    //     for (int i=0; i<2; i++){
    //         aqIVP1.set_time_steps(steps[i]);
    //         aqIVP1.reset_f_evaluations();
    //         switch (j) {
    //             case 0:
    //                 aqIVP1.solve_rungekutta(4);
    //                 break;
    //             case 1:
    //                 aqIVP1.solve_rungekutta(5);
    //                 break;
    //             default:
    //                 printf("Undifined!\n");
    //                 return;
    //         }
    //         evaluations = aqIVP1.get_f_evaluations();
    //         aqIVP1.solve_integral();
    //         aqIVP1.calculate_exact_error();
    //         error_list = aqIVP1.get_y_cumulative_error();
    //         error_end = error_list[steps[i]];
    //         error_max = max_array(error_list, steps[i]+1);
    //         // printf("%10d\t%16d\t%16.10g\t%16.10g\n",steps[i], evaluations, error_end, error_max);
    //         fprintf(outFile,"%10d\t%16d\t%16.10g\t%16.10g\n",steps[i], evaluations, error_end, error_max);
    //     }
    //     fclose(outFile);
    // }
}

int main(){
    #ifdef _WIN32
        system("cls");
    #elif __unix__
        system("clear");
    #else
        std::cout << "Running on an unknown system" << std::endl;
    #endif
    tests_rungekutta();
    Fetkovich_tests();
}