#include <iostream>
#include <cmath>
#include <limits>
#include <cstdio>
#include <string>
#include <cstdlib>
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
class Spline{
    public:
        Spline();
        void add_points(double *x, double *y, int number_points);
        void add_y_prime(double y_prime_initial, double y_prime_end);
        void build_natural();
        void build_fixed();

        void set_extrapolation(bool enable_extrapolation);
        double get_y(double x) const;
        double get_y_prime(double x) const;
        double get_y_prime2(double x) const;
        double get_int_y(double x) const;
        double get_int_y(double x_initial, double x_end) const;

    private:
        bool has_data() const;
        bool has_spline() const;
        bool has_primes() const;
        void reset_data();
        void reset_spline();

        bool is_extrapolation(double x) const;
        int find_interval(double x) const;
        double get_int_y_point(double x, int interval) const;

        int n_points;
        double *x_pointer;
        double *a_pointer;
        double *b_pointer;
        double *c_pointer;
        double *d_pointer;
        double *y_prime_pointer;
        bool can_extrapolate;
};

Spline::Spline():
    n_points(0),
    x_pointer(nullptr), a_pointer(nullptr),
    b_pointer(nullptr), c_pointer(nullptr),
    d_pointer(nullptr), y_prime_pointer(nullptr),
    can_extrapolate(false)
    {};

void Spline::reset_data(){
    n_points = 0;
    x_pointer = nullptr;
    a_pointer = nullptr;
    y_prime_pointer = nullptr;
}

void Spline::reset_spline(){
    b_pointer = nullptr;
    c_pointer = nullptr;
    d_pointer = nullptr;
}

void Spline::add_points(double *x, double *y, int number_points){
    reset_spline();
    if (number_points<3){
        printf("At least 3 points are needed to build a spline.\n");
        return;
    }
    for(int i=0; i<number_points-1; i++){
        if (x[i+1] <= x[i]){
            printf("X values must be in ascending order. Cannot continue.\n");
            return;
        };
    }
    x_pointer = copy_array(x, number_points);
    a_pointer = copy_array(y, number_points);
    n_points = number_points;
    y_prime_pointer = nullptr;
}

void Spline::add_y_prime(double y_prime_initial, double y_prime_end){
    if (!has_data()){
        printf("No data found. Inform (x,y) values first.\n");
        return;
    }
    double* y_prime = new double[n_points];
    for (int i=1; i<n_points-1;i++){
        y_prime[i] = 0;
    }
    y_prime[0] = y_prime_initial;
    y_prime[n_points-1] = y_prime_end;
    y_prime_pointer = y_prime;
}

bool Spline::has_data() const{
    return x_pointer != nullptr && a_pointer != nullptr && n_points > 2;
}

bool Spline::has_spline() const{
    return x_pointer != nullptr && a_pointer != nullptr && b_pointer != nullptr && c_pointer != nullptr && d_pointer != nullptr;
}

bool Spline::has_primes() const{
    return y_prime_pointer != nullptr;
}

void Spline::build_natural(){
    reset_spline();
    if (!has_data()){
        printf("No data found. Cannot build natural spline.\n");
        return;
    }
    double* b = new double[n_points];
    double* c = new double[n_points];
    double* d = new double[n_points];

    double* h = new double[n_points];
    double* r = new double[n_points];
    double* m = new double[n_points];
    double* z = new double[n_points];
    double* alpha = new double[n_points];

    int n = n_points - 1;
    for(int i=0; i<n; i++){
        h[i] = x_pointer[i+1] - x_pointer[i];
    }

    for(int i=1; i<n; i++){
        alpha[i] = 3. / h[i] * (a_pointer[i+1] - a_pointer[i]) - 3. / h[i-1] * (a_pointer[i] - a_pointer[i-1]);
    }
    r[0] = 1.;
    m[0] = 0.;
    z[0] = 0.;

    for(int i=1; i<n; i++){
        r[i] = 2. * (x_pointer[i+1] - x_pointer[i-1]) - h[i-1] * m[i-1];
        if (fabs(r[i]) < eps){
            printf("Error in natural spline calculation: |r[i]| too small. Check data around point #%i. Cannot continue.\n",i);
            return;
        }
        m[i] = h[i] / r[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / r[i];
    }

    r[n] = 1.;
    z[n] = 0.;
    c[n] = 0.;

    for(int i=n-1; i>-1; i--){
        c[i] = z[i] - m[i] * c[i+1];
        b[i] = (a_pointer[i+1] - a_pointer[i]) / h[i] - h[i] * (c[i+1] + 2*c[i]) / 3;
        d[i] = (c[i+1] - c[i]) / 3 / h[i];
    }

    if (fabs(c[0]) > eps){
        printf("Natural spline calculation failed: c[0] != 0. Check code. Cannot continue.\n");
        return;
    };

    b_pointer = b;
    c_pointer = c;
    d_pointer = d;

    if (fabs(get_y_prime2(x_pointer[0])) > eps){
        printf("Natural spline calculation failed: y''[0] != 0. Check code. Cannot continue.\n");
        reset_spline();
        return;
    }
    if (fabs(get_y_prime2(x_pointer[n_points-1])) > eps){
        printf("Natural spline calculation failed: y''[n] != 0. Check code. Cannot continue.\n");
        reset_spline();
        return;
    }
}

void Spline::build_fixed(){
    reset_spline();
    if (!has_data()){
        printf("No data found. Cannot build fixed spline.\n");
        return;
    }
    if (!has_primes()){
        printf("No y prime data found. Cannot build fixed spline.\n");
        return;
    }
    double* b = new double[n_points];
    double* c = new double[n_points];
    double* d = new double[n_points];

    double* h = new double[n_points];
    double* r = new double[n_points];
    double* m = new double[n_points];
    double* z = new double[n_points];
    double* alpha = new double[n_points];

    int n = n_points - 1;
    for(int i=0; i<n; i++){
        h[i] = x_pointer[i+1] - x_pointer[i];
    }

    alpha[0] = 3. / h[0] * (a_pointer[1] - a_pointer[0]) - 3. * y_prime_pointer[0];
    alpha[n] = 3. * y_prime_pointer[n] - 3. / h[n-1] * (a_pointer[n] - a_pointer[n-1]);
    for(int i=1; i<n; i++){
        alpha[i] = 3. / h[i] * (a_pointer[i+1] - a_pointer[i]) - 3. / h[i-1] * (a_pointer[i] - a_pointer[i-1]);
    }
    r[0] = 2*h[0];
    m[0] = 0.5;
    z[0] = alpha[0] / r[0];

    for(int i=1; i<n; i++){
        r[i] = 2. * (x_pointer[i+1] - x_pointer[i-1]) - h[i-1] * m[i-1];
        if (fabs(r[i]) < eps){
            printf("Error in fixed spline calculation: |r[i]| too small. Check data around point #%i. Cannot continue.\n",i);
            return;
        }
        m[i] = h[i] / r[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / r[i];
    }

    r[n] = h[n-1] * (2 - m[n-1]);
    z[n] = (alpha[n] - h[n-1] * z[n-1]) / r[n];
    c[n] = z[n];

    for(int i=n-1; i>-1; i--){
        c[i] = z[i] - m[i] * c[i+1];
        b[i] = (a_pointer[i+1] - a_pointer[i]) / h[i] - h[i] * (c[i+1] + 2*c[i]) / 3;
        d[i] = (c[i+1] - c[i]) / 3 / h[i];
    }

    if (fabs(b[0] - y_prime_pointer[0]) > eps){
        printf("Fixed spline calculation failed: b[0] != y_prime_initial. Check code. Cannot continue.\n");
        return;
    };

    b_pointer = b;
    c_pointer = c;
    d_pointer = d;

    if (fabs(get_y_prime(x_pointer[n_points-1]) - y_prime_pointer[n]) > eps){
        printf("Fixed spline calculation failed: y'[n] != y_prime_end. Check code. Cannot continue.\n");
        reset_spline();
        return;
    }
}

void Spline::set_extrapolation(bool enable_extrapolation){
    can_extrapolate = enable_extrapolation;
}
bool Spline::is_extrapolation(double x) const{
    if (!has_data()){
        return false;
    }
    return (x < x_pointer[0]) || (x > x_pointer[n_points-1]);
}

int Spline::find_interval(double x) const{
    if (!has_data()){
        printf("No data found.\n");
        return -1;
    }
    for (int i=1; i<n_points-1; i++){
        if( x <= x_pointer[i]){
            return i-1;
        }
    }
    return n_points-2;
}

double Spline::get_y(double x) const{
    if (!has_spline()){
        printf("Spline not built.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (!can_extrapolate && is_extrapolation(x)){
        return std::numeric_limits<double>::quiet_NaN();
    }
    int j = find_interval(x);
    double h = x - x_pointer[j];

    double y = a_pointer[j];
    double dx = h;
    y += b_pointer[j] * dx;
    dx *= h;
    y += c_pointer[j] * dx;
    dx *= h;
    y += d_pointer[j] * dx;

    return y;
}
double Spline::get_y_prime(double x) const{
    if (!has_spline()){
        printf("Spline not built.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (!can_extrapolate && is_extrapolation(x)){
        return std::numeric_limits<double>::quiet_NaN();
    }
    int j = find_interval(x);
    double h = x - x_pointer[j];

    double y = b_pointer[j];
    double dx = h;
    y += c_pointer[j] * dx * 2.;
    dx *= h;
    y += d_pointer[j] * dx * 3.;

    return y;
}
double Spline::get_y_prime2(double x) const{
    if (!has_spline()){
        printf("Spline not built.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (!can_extrapolate && is_extrapolation(x)){
        return std::numeric_limits<double>::quiet_NaN();
    }
    int j = find_interval(x);
    double h = x - x_pointer[j];

    double y = c_pointer[j] * 2.;
    y += d_pointer[j] * h * 6.;

    return y;
}

double Spline::get_int_y_point(double x, int interval) const{
    int j = interval;
    double h = x - x_pointer[j];

    double dx = h;
    double y = a_pointer[j] * dx;
    dx *= h;
    y += b_pointer[j] * dx / 2.;
    dx *= h;
    y += c_pointer[j] * dx / 3.;
    dx *= h;
    y += d_pointer[j] * dx / 4.;

    return y;
}
double Spline::get_int_y(double x) const{
    if (!has_data()){
        printf("No data found.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return get_int_y(x_pointer[0], x);
}
double Spline::get_int_y(double x_initial, double x_end) const{
    if (!has_spline()){
        printf("Spline not built.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (!can_extrapolate && (is_extrapolation(x_initial) || is_extrapolation(x_end))){
        return std::numeric_limits<double>::quiet_NaN();
    }
    int j_initial;
    int j_end;
    double y;

    j_initial = find_interval(x_initial);
    j_end = find_interval(x_end);

    y = -get_int_y_point(x_initial, j_initial);

    if( j1 != j0){
        for( int j=j_initial; j<j_end; j++){
            y += get_int_y_point(x_pointer[j+1], j);
        }
    }
    y += get_int_y_point(x_end, j_end);

    return y;
}

class Integration{
    public:
        Integration();
        void add_points(double *x, double *y, int number_points);
        double solve_Simpson() const;
        double solve_Trapezoidal() const;
        double* solve() const;

    private:
        bool has_data() const;

        int n_points;
        double *x_pointer;
        double *y_pointer;
};

Integration::Integration():
    n_points(0),
    x_pointer(nullptr), y_pointer(nullptr)
    {};

void Integration::add_points(double *x, double *y, int number_points){
    if (number_points<2){
        printf("At least 2 points are needed to compute an integral.\n");
        return;
    }
    for(int i=0; i<number_points-1; i++){
        if (x[i+1] <= x[i]){
            printf("X values must be in ascending order. Cannot continue.\n");
            return;
        };
    }
    x_pointer = copy_array(x, number_points);
    y_pointer = copy_array(y, number_points);
    n_points = number_points;
}

bool Integration::has_data() const{
    return x_pointer != nullptr && y_pointer != nullptr && n_points > 2;
}

double Integration::solve_Simpson() const{
    if (!has_data()){
        printf("No data defined yet.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (n_points % 2 == 0){
        printf("Simpson only works with an odd number of points.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    double h = (x_pointer[n_points-1] - x_pointer[0]) / (n_points-1);
    for (int i=1; i<n_points; i++){
        if (abs(h - (x_pointer[i] - x_pointer[i-1])) > eps){
            printf("Simpson only works with a constant step.\n");
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    double s0 = y_pointer[n_points-1] + y_pointer[0];
    double s1 = 0.;
    double s2 = 0.;
    for (int i=1; i<n_points-1; i++){
        if (i % 2 == 0){
            s2 += y_pointer[i];
        } else{
            s1 += y_pointer[i];
        }
    }
    return h/3. * (s0 + 2.*s2 + 4.*s1);
}

double Integration::solve_Trapezoidal() const{
    if (!has_data()){
        printf("No data defined yet.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    double h = (x_pointer[n_points-1] - x_pointer[0]) / (n_points-1);
    for (int i=1; i<n_points; i++){
        if (abs(h - (x_pointer[i] - x_pointer[i-1])) > eps){
            printf("Trapezoidal rule integration only works with a constant step.\n");
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    double s0 = y_pointer[n_points-1] + y_pointer[0];
    double s1 = 0.;
    for (int i=1; i<n_points-1; i++){
        s1 += y_pointer[i];
    }
    return h/2. * (s0 + 2.*s1);
}

double* Integration::solve() const{
    if (!has_data()){
        printf("No data defined yet.\n");
        return nan_array(1);
    }
    double* h = new double[n_points-1];
    for (int i=1; i<n_points; i++){
        h[i-1] = x_pointer[i] - x_pointer[i-1];
        if (h[i-1] <= 0){
            printf("X-values must be in ascending order.\n");
            return nan_array(n_points);
        }
    }
    double* s = new double[n_points];
    s[0] = 0.;
    double w1, w2, w3;
    for (int i=1; i<n_points; i++){
        if (i % 2 == 1){
            s[i] = s[i-1] + h[i-1]/2. * (y_pointer[i] + y_pointer[i-1]);
        } else{
            w1 = 2.*h[i-2] + h[i-1] * (1. - h[i-1]/h[i-2]);
            w2 = pow(h[i-2] + h[i-1], 3.) / (h[i-2]*h[i-1]);
            w3 = 2.*h[i-1] + h[i-2] * (1. - h[i-2]/h[i-1]);
            s[i] = s[i-2] + 1./6. * (y_pointer[i-2]*w1 + y_pointer[i-1]*w2 + y_pointer[i]*w3);
        }
    }
    return s;
}

void tests_integration(){
    Integration integration;
    int n = 18;
    double* x = new double[n+1];
    double* y = new double[n+1];
    for (int i=0; i<n+1; i++){
        x[i] = M_PI * i/n;
        y[i] = sin(x[i]);
    }
    integration.add_points(x,y,n+1);
    printf("Integration of Sin(x) from 0 to pi\n");
    printf("  Simpson: %g (expected 2.0000104)\n", integration.solve_Simpson());
    printf("  Trapezoidal rule: %g (expected 1.9949205)\n", integration.solve_Trapezoidal());
    double* s = integration.solve();
    printf("  'General' integrator:\n");
    printf("%5s\t%10s\t%10s\t%10s\t%10s\t\n","i","x","Aprox","Exact","Error");
    for (int i=0; i<n+1; i++){
        printf("%5d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t\n",i,x[i],s[i],1-cos(x[i]),abs(s[i]-1+cos(x[i])));
    }
}


class IVP{
    public:
        IVP();
        void set_f(double (*f)(double, double));
        void set_exact(double (*f)(double));
        void set_exact_cumulative(double (*f)(double));
        void set_dfdt(double (*f)(double, double));
        void set_delfdely(double (*f)(double, double));
        void set_y_initial(double y);
        void set_t_initial(double t);
        void set_t_end(double t);
        void set_time_steps(int n);
        void set_relative_error(bool relative);
        void set_extra_search_points(int n);

        void set_Lipschitz_L(double L);
        void estimate_Lipschitz();

        void set_max_dfdt(double M);
        void estimate_max_dfdt();

        double get_f(double t, double y);
        double get_dfdt(double t, double y) const;
        double get_delfdely(double t, double y) const;
        double get_exact(double t) const;
        double get_exact_cumulative(double t) const;
        double get_Lipschitz() const;
        double get_max_dfdt() const;

        double* get_t() const;
        double* get_y() const;
        double* get_y_cumulative() const;
        double* get_y_error() const;
        double* get_y_cumulative_error() const;
        int get_f_evaluations() const;

        void reset_f_evaluations();
        void calculate_exact_error();
        void print_solution();
        void print_solution(string filename);

        void initialize_problem();

        void solve_euler();
        void estimate_error_euler();
        void solve_euler_aitken();
        void solve_euler_aitken(int n1, int n2);

        void solve_rungekutta();
        void solve_rungekutta(int n);
        void solve_rungekutta_aitken();
        void solve_rungekutta_aitken(int n);
        void solve_rungekutta_aitken(int n, int n1, int n2);

        void solve_adams();
        void solve_adams(int steps, bool implicit);
        void set_adams_convergence_limit(double value);
        void set_adams_convergence_iterartions(int value);
        void set_adams_convergence_print(bool value);

        void solve_integral();

    private:
        void solve_rungekutta2();
        void solve_rungekutta3();
        void solve_rungekutta4();
        void solve_rungekutta5();
        void solve_rungekutta2(bool one_step, int current_step);
        void solve_rungekutta3(bool one_step, int current_step);
        void solve_rungekutta4(bool one_step, int current_step);
        void solve_rungekutta5(bool one_step, int current_step);

        void solve_adams_exp_2(bool one_step, int current_step);
        void solve_adams_exp_3(bool one_step, int current_step);
        void solve_adams_exp_4(bool one_step, int current_step);
        void solve_adams_exp_5(bool one_step, int current_step);
        void solve_adams_1(bool one_step, int current_step);
        void solve_adams_2(bool one_step, int current_step);
        void solve_adams_3(bool one_step, int current_step);
        void solve_adams_4(bool one_step, int current_step);

        void reset_results();
        void reset_error_estimate();

        bool has_f() const;
        bool has_exact() const;
        bool has_exact_cumulative() const;
        bool has_solution() const;
        bool has_solution_cumulative() const;
        bool has_exact_error() const;
        bool has_estimated_error() const;
        bool has_dfdt() const;
        bool has_delfdely() const;
        bool has_Lipschitz() const;
        bool has_max_dfdt() const;

        double (*f_pointer)(double, double);
        double (*exact_pointer)(double);
        double (*exact_cum_pointer)(double);
        double (*dfdt_pointer)(double, double);
        double (*delfdely_pointer)(double, double);
        double y_initial;
        double t_initial;
        double t_end;
        int time_steps;
        double h_current;
        double *t_pointer;
        double *y_pointer;
        double *y_exact_pointer;
        double *y_error_pointer;
        double *y_error_limit_pointer;
        double *y_cum_pointer;
        double *y_exact_cum_pointer;
        double *y_error_cum_pointer;
        int *evaluations_pointer;
        int extra_search_points;
        bool relative_error;
        double Lipschitz;
        bool calculated_L;
        double max_dfdt;
        bool calculated_M;
        int f_evaluations;
        double conv_adams;
        int max_iter_adams;
        bool print_conv_adams;
};

IVP::IVP():
    f_pointer(nullptr), exact_pointer(nullptr),
    exact_cum_pointer(nullptr),
    dfdt_pointer(nullptr),
    delfdely_pointer(nullptr),
    y_initial(0), t_initial(0),
    t_end(1), time_steps(100), h_current(1.),
    t_pointer(nullptr), y_pointer(nullptr),
    y_exact_pointer(nullptr), y_error_pointer(nullptr),
    y_error_limit_pointer(nullptr),
    y_cum_pointer(nullptr),
    y_exact_cum_pointer(nullptr),
    y_error_cum_pointer(nullptr),
    evaluations_pointer(nullptr),
    extra_search_points(20), relative_error(false),
    Lipschitz(0), calculated_L(false),
    max_dfdt(0), calculated_M(false),
    f_evaluations(0),
    conv_adams(1e-3), max_iter_adams(1), print_conv_adams(false)
    {};

void IVP::reset_results(){
    t_pointer = nullptr;
    y_pointer = nullptr;
    y_cum_pointer = nullptr;
    calculated_L = false;
    calculated_M = false;
    // reset_f_evaluations();
    reset_error_estimate();
}
void IVP::reset_f_evaluations(){
    f_evaluations = 0;
}
void IVP::reset_error_estimate(){
    y_exact_pointer = nullptr;
    y_error_pointer = nullptr;
    y_error_cum_pointer = nullptr;
    y_error_limit_pointer = nullptr;
}
void IVP::set_f(double (*f)(double, double)){
    f_pointer = f;
    reset_results();
}
void IVP::set_dfdt(double (*f)(double, double)){
    dfdt_pointer = f;
    calculated_M = false;
}
void IVP::set_delfdely(double (*f)(double, double)){
    delfdely_pointer = f;
    calculated_L = false;
}
void IVP::set_exact(double (*f)(double)){
    exact_pointer = f;
}
void IVP::set_exact_cumulative(double (*f)(double)){
    exact_cum_pointer = f;
}
void IVP::set_y_initial(double y){
    y_initial = y;
    reset_results();
}
void IVP::set_t_initial(double t){
    t_initial = t;
    if (abs(t-t_end)<eps){
        t_end = t_initial + 1;
        printf("Cannot set initial time equal to end time. Setting end time to initial time + 1: %g\n", t_end);
    }
    reset_results();
}
void IVP::set_t_end(double t){
    if (abs(t-t_initial)<eps){
        printf("Cannot set end time equal to initial time. Keeping previous value: %g\n", t_end);
    } else{
        t_end = t;
    }
    reset_results();
}
void IVP::set_time_steps(int n){
    if (n<1){
        printf("Cannot set less than one time-step. Keeping previous value: %i\n", time_steps);
    } else{
        time_steps = n;
    }
    reset_results();
}
void IVP::set_relative_error(bool relative){
    if (relative_error != relative){
        relative_error = relative;
        y_error_pointer = nullptr;
        y_error_cum_pointer = nullptr;
    }
}
void IVP::set_Lipschitz_L(double L){
    if (L<=0){
        printf("Cannot set a non-positive value.\n");
        calculated_L = false;
    } else{
        Lipschitz = L;
        calculated_L = true;
    }
}
void IVP::set_max_dfdt(double M){
    if (M<=0){
        printf("Cannot set a non-positive value.\n");
        calculated_M = false;
    } else{
        max_dfdt = M;
        calculated_M = true;
    }
}

void IVP::set_extra_search_points(int n){
    if (n<0){
        printf("Cannot set a negative value. Assigning zero.");
        extra_search_points = 0;
    } else{
        extra_search_points = n;
    }
}

double IVP::get_f(double t, double y){
    if (!has_f()){
        printf("Function f not defined.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    f_evaluations++;
    return f_pointer(t,y);
}
double IVP::get_dfdt(double t, double y) const{
    if (!has_dfdt()){
        printf("Function df/dt not defined.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return dfdt_pointer(t,y);
}
double IVP::get_delfdely(double t, double y) const{
    if (!has_delfdely()){
        printf("Function del_f/del_y not defined.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return delfdely_pointer(t,y);
}
double IVP::get_exact(double t) const{
    if (!has_exact()){
        printf("Exact response not defined.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return exact_pointer(t);
}
double IVP::get_exact_cumulative(double t) const{
    if (!has_exact_cumulative()){
        printf("Exact cumulative response not defined.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return exact_cum_pointer(t);
}
double IVP::get_Lipschitz() const{
    if (!calculated_L){
        printf("Lipschitz L was neither provided nor estimated.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return Lipschitz;
}
void IVP::estimate_Lipschitz(){
    if (!has_solution()){
        printf("No solutions found. Cannot continue.\n");
        return;
    }
    if (!has_delfdely()){
        printf("Function del_f/del_y is not defined. Cannot continue.\n");
        return;
    }

    double L=0;
    double Lmax=abs(get_delfdely(t_pointer[0], y_pointer[0]));
    double y;
    double t;

    Spline spl;
    if (extra_search_points > 0){
        spl.add_points(t_pointer,y_pointer,time_steps+1);
        if (has_dfdt()){
            spl.add_y_prime(get_dfdt(t_pointer[0], y_pointer[0]), get_dfdt(t_pointer[time_steps], y_pointer[time_steps]));
            spl.build_fixed();
        } else{
            spl.build_natural();
        }
    }

    for (int i=1; i<time_steps+1; i++){
        t = t_pointer[i];
        y = y_pointer[i];
        L = abs(get_delfdely(t, y));
        if (L>Lmax){
            Lmax = L;
        }

        for (int j=0; j<extra_search_points; j++){
            t = t_pointer[i-1] + (t_pointer[i] - t_pointer[i-1]) * (j+1) / (extra_search_points+1);
            y = spl.get_y(t);
            L = abs(get_delfdely(t, y));
            if (L>Lmax){
                Lmax = L;
            }
        }
    }
    Lipschitz = Lmax;
    calculated_L = true;
}

double IVP::get_max_dfdt() const{
    if (!calculated_M){
        printf("Max(df/dt) was neither provided nor estimated.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return max_dfdt;
}
void IVP::estimate_max_dfdt(){
    if (!has_solution()){
        printf("No solutions found. Cannot continue.\n");
        return;
    }
    if (!has_dfdt()){
        printf("Function df/dt is not defined. Cannot continue.\n");
        return;
    }

    double M=0;
    double Mmax=abs(get_dfdt(t_pointer[0], y_pointer[0]));
    double y;
    double t;

    Spline spl;
    if (extra_search_points > 0){
        spl.add_points(t_pointer,y_pointer,time_steps+1);
        spl.add_y_prime(get_dfdt(t_pointer[0], y_pointer[0]), get_dfdt(t_pointer[time_steps], y_pointer[time_steps]));
        spl.build_fixed();
    }

    for (int i=1; i<time_steps+1; i++){
        t = t_pointer[i];
        y = y_pointer[i];
        M = abs(get_dfdt(t, y));
        if (M>Mmax){
            Mmax = M;
        }

        for (int j=0; j<extra_search_points; j++){
            t = t_pointer[i-1] + (t_pointer[i] - t_pointer[i-1]) * (j+1) / (extra_search_points+1);
            y = spl.get_y(t);
            M = abs(get_dfdt(t, y));
            if (M>Mmax){
                Mmax = M;
            }
        }
    }
    max_dfdt = Mmax;
    calculated_M = true;
}

bool IVP::has_f() const{
    return f_pointer != nullptr;
}
bool IVP::has_exact() const{
    return exact_pointer != nullptr;
}
bool IVP::has_exact_cumulative() const{
    return exact_cum_pointer != nullptr;
}
bool IVP::has_solution() const{
    return t_pointer != nullptr && y_pointer != nullptr;
}
bool IVP::has_solution_cumulative() const{
    return y_cum_pointer != nullptr;
}
bool IVP::has_exact_error() const{
    return y_exact_pointer != nullptr && y_error_pointer != nullptr;
}
bool IVP::has_estimated_error() const{
    return y_error_limit_pointer != nullptr;
}
bool IVP::has_dfdt() const{
    return dfdt_pointer != nullptr;
}
bool IVP::has_delfdely() const{
    return delfdely_pointer != nullptr;
}
bool IVP::has_Lipschitz() const{
    return calculated_L;
}
bool IVP::has_max_dfdt() const{
    return calculated_M;
}

double* IVP::get_t() const{
    return copy_array(t_pointer, time_steps+1);
}
double* IVP::get_y() const{
    return copy_array(y_pointer, time_steps+1);
}
double* IVP::get_y_cumulative() const{
    return copy_array(y_cum_pointer, time_steps+1);
}
double* IVP::get_y_error() const{
    return copy_array(y_error_pointer, time_steps+1);
}
double* IVP::get_y_cumulative_error() const{
    return copy_array(y_error_cum_pointer, time_steps+1);
}

int IVP::get_f_evaluations() const{
    return f_evaluations;
}

void IVP::calculate_exact_error(){
    if (!has_solution()){
        printf("No solutions found. Cannot continue.\n");
        return;
    }
    if (!has_exact()){
        printf("Exact solution not defined. Cannot continue.\n");
        return;
    }
    double* y_exact = new double[time_steps+1];
    double* y_error = new double[time_steps+1];
    for (int i=0; i<time_steps+1; i++){
        y_exact[i] = get_exact(t_pointer[i]);
        y_error[i] = abs(y_exact[i] - y_pointer[i]);
        if (relative_error){
            if (abs(y_exact[i]) > eps){
                y_error[i] = y_error[i] / abs(y_exact[i]);
            } else{
                y_error[i] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
    y_exact_pointer = y_exact;
    y_error_pointer = y_error;

    if (has_exact_cumulative() && has_solution_cumulative()){
        double* y_exact_cum = new double[time_steps+1];
        double* y_error_cum = new double[time_steps+1];
        for (int i=0; i<time_steps+1; i++){
            y_exact_cum[i] = get_exact_cumulative(t_pointer[i]);
            y_error_cum[i] = abs(y_exact_cum[i] - y_cum_pointer[i]);
            if (relative_error){
                if (abs(y_exact_cum[i]) > eps){
                    y_error_cum[i] = y_error_cum[i] / abs(y_exact_cum[i]);
                } else{
                    y_error_cum[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
        y_exact_cum_pointer = y_exact_cum;
        y_error_cum_pointer = y_error_cum;
    }
}

void IVP::print_solution(){
    if (!has_solution()){
        printf("There is no solution to print.\n");
        return;
    }

    printf("%5s\t%16s\t%16s","i","t","y_aprox");
    if (has_estimated_error()){
        printf("\t%16s","error_limit");
    }
    if (has_exact_error()){
        printf("\t%16s\t%16s","y_exact","y_error");
    }
    printf("\n");

    for (int i=0; i<time_steps+1; i++){
        printf("%5i\t%16.10g\t%16.10g", i, t_pointer[i], y_pointer[i]);
        if (has_estimated_error()){
            printf("\t%16.10g", y_error_limit_pointer[i]);
        }
        if (has_exact_error()){
            printf("\t%16.10g\t%16.10g", y_exact_pointer[i],y_error_pointer[i]);
        }
        printf("\n");
    }
    if (relative_error){
        printf("* Errors are relative.\n");
    }
}

void IVP::print_solution(string filename){
    if (!has_solution()){
        printf("There is no solution to print.\n");
        return;
    }

    FILE* outFile = fopen(filename.c_str(), "w");
    if (outFile == nullptr) {
        printf("Coud not create file: %s\n",filename.c_str());
        return;
    }

    fprintf(outFile,"%5s\t%16s\t%16s\t%16s","i","t","evaluations","y_aprox");
    if (has_estimated_error()){
        fprintf(outFile,"\t%16s","error_limit");
    }
    if (has_exact_error()){
        fprintf(outFile,"\t%16s\t%16s","y_exact","exact_error");
    }
    if (has_solution_cumulative()){
        fprintf(outFile,"\t%16s","Sy_aprox");
        if (has_exact_cumulative()){
            fprintf(outFile,"\t%16s\t%16s","Sy_exact","S_error");
        }
    }
    fprintf(outFile,"\n");

    for (int i=0; i<time_steps+1; i++){
        fprintf(outFile,"%5i\t%16.10g\t%16i\t%16.10g", i, t_pointer[i], evaluations_pointer[i], y_pointer[i]);
        if (has_estimated_error()){
            fprintf(outFile,"\t%16.10g", y_error_limit_pointer[i]);
        }
        if (has_exact_error()){
            fprintf(outFile,"\t%16.10g\t%16.10g", y_exact_pointer[i],y_error_pointer[i]);
        }
        if (has_solution_cumulative()){
            fprintf(outFile,"\t%16.10g",y_cum_pointer[i]);
            if (has_exact_cumulative()){
                fprintf(outFile,"\t%16.10g\t%16.10g",y_exact_cum_pointer[i],y_error_cum_pointer[i]);
            }
        }
        fprintf(outFile,"\n");
    }
    fclose(outFile);
}

void IVP::initialize_problem(){
    double* t = new double[time_steps+1];
    double* y = new double[time_steps+1];
    int* f = new int[time_steps+1];
    t[0] = t_initial;
    y[0] = y_initial;
    f[0] = 0;
    t_pointer = t;
    y_pointer = y;
    evaluations_pointer = f;
    h_current = (t_end - t_initial) / time_steps;
}

void IVP::solve_euler(){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    initialize_problem();
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    for (int i=0; i<time_steps; i++){
        y[i+1] = y[i] + h * f_pointer(t[i], y[i]);
        t[i+1] = t[0] + (i+1)*h;
        f_evaluations++;
        f[i+1] = f_evaluations;
    }
    reset_error_estimate();
}

void IVP::estimate_error_euler(){
    if (!has_solution()){
        printf("No solution found. Cannot continue.\n");
        return;
    }
    if (!has_Lipschitz()){
        printf("Lipschitz L neither provided nor estimated. Cannot continue.\n");
        return;
    }
    if (!has_max_dfdt()){
        printf("Max(df/dt) neither provided nor estimated. Cannot continue.\n");
        return;
    }

    double L = get_Lipschitz();
    double M = get_max_dfdt();
    double h = (t_end - t_initial) / time_steps;

    double* error_limit = new double[time_steps+1];

    error_limit[0] = 0;

    for (int i=1; i<time_steps+1; i++){
        error_limit[i] = h*M/(2*L)*(exp(L*(t_pointer[i] - t_pointer[0])) - 1);
        if (relative_error){
            if (abs(y_pointer[i]) > eps){
                error_limit[i] = error_limit[i] / abs(y_pointer[i]);
            } else{
                error_limit[i] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }

    y_error_limit_pointer = error_limit;
}

void IVP::solve_euler_aitken(){
    solve_euler_aitken(2,4);
}

void IVP::solve_euler_aitken(int n1, int n2){
    int time_steps_original = time_steps;

    set_time_steps(n2 * time_steps_original);
    solve_euler();
    double* y2 = copy_array(y_pointer, time_steps+1);

    set_time_steps(n1 * time_steps_original);
    solve_euler();
    double* y1 = copy_array(y_pointer, time_steps+1);

    set_time_steps(time_steps_original);
    solve_euler();
    double* y0 = copy_array(y_pointer, time_steps+1);

    for (int i=1; i<time_steps+1; i++){
        y_pointer[i] = y0[i] - pow(y1[n1*i] - y0[i], 2) / (y2[n2*i] - 2*y1[n1*i] + y0[i]);
    }
}

void IVP::solve_rungekutta(){
    solve_rungekutta4();
}

void IVP::solve_rungekutta(int n){
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
        printf("No Runge-Kutta methods higher than order 5 were implemented.\n");
        break;
    }
}

void IVP::solve_rungekutta2(){
    initialize_problem();
    solve_rungekutta2(false, 1);
}

void IVP::solve_rungekutta2(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }
    double k1,k2;
    for (int i=(current_step-1); i<steps; i++){
        t[i+1] = t[0] + (i+1)*h;
        k1 = h * f_pointer(t[i], y[i]);
        k2 = h * f_pointer(t[i]+h/2., y[i]+k1/2.);
        f_evaluations += 2;
        f[i+1] = f_evaluations;

        y[i+1] = y[i] + k2;
    }
    reset_error_estimate();
}

void IVP::solve_rungekutta3(){
    initialize_problem();
    solve_rungekutta3(false, 1);
}

void IVP::solve_rungekutta3(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }
    double k1,k2,k3;
    for (int i=(current_step-1); i<steps; i++){
        t[i+1] = t[0] + (i+1)*h;
        k1 = h * f_pointer(t[i], y[i]);
        k2 = h * f_pointer(t[i]+h/3., y[i]+k1/3.);
        k3 = h * f_pointer(t[i]+h*2./3., y[i]+k2*2./3.);
        f_evaluations += 3;
        f[i+1] = f_evaluations;

        y[i+1] = y[i] + 1./4. * (k1 + 3.*k3);
    }
    reset_error_estimate();
}

void IVP::solve_rungekutta4(){
    initialize_problem();
    solve_rungekutta4(false, 1);
}

void IVP::solve_rungekutta4(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }
    double k1,k2,k3,k4;
    for (int i=(current_step-1); i<steps; i++){
        t[i+1] = t[0] + (i+1)*h;
        k1 = h * f_pointer(t[i], y[i]);
        k2 = h * f_pointer(t[i]+h/2., y[i]+k1/2.);
        k3 = h * f_pointer(t[i]+h/2., y[i]+k2/2.);
        k4 = h * f_pointer(t[i+1], y[i]+k3);
        f_evaluations += 4;
        f[i+1] = f_evaluations;

        y[i+1] = y[i] + 1./6. * (k1 + 2.*k2 + 2.*k3 + k4);
    }
    reset_error_estimate();
}

void IVP::solve_rungekutta5(){
    initialize_problem();
    solve_rungekutta5(false, 1);
}

void IVP::solve_rungekutta5(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }
    double k1,k2,k3,k4,k5,k6;
    for (int i=(current_step-1); i<steps; i++){
        t[i+1] = t[0] + (i+1)*h;
        k1 = h * f_pointer(t[i]          , y[i]);
        k2 = h * f_pointer(t[i]+h/4.     , y[i] + k1/4.);
        k3 = h * f_pointer(t[i]+h*3./8.  , y[i] + k1*3./32.      + k2*9./32.);
        k4 = h * f_pointer(t[i]+h*12./13., y[i] + k1*1932./2197. - k2*7200./2197. + k3*7296./2197.);
        k5 = h * f_pointer(t[i]+h        , y[i] + k1*439./216.   - k2*8.          + k3*3680./513.  - k4*845./4104.);
        k6 = h * f_pointer(t[i]+h/2.     , y[i] - k1*8./27.      + k2*2.          - k3*3544./2565. + k4*1859./4104. - k5*11./40.);
        f_evaluations += 6;
        f[i+1] = f_evaluations;

        y[i+1] = y[i] + 16./135.*k1 + 6656./12825.*k3 + 28561./56430.*k4 - 9./50.*k5 + 2./55.*k6;
    }
    reset_error_estimate();
}

void IVP::solve_rungekutta_aitken(){
    solve_rungekutta_aitken(4, 2, 4);
}
void IVP::solve_rungekutta_aitken(int n){
    solve_rungekutta_aitken(n, 2, 4);
}
void IVP::solve_rungekutta_aitken(int n, int n1, int n2){
    int time_steps_original = time_steps;

    set_time_steps(n2 * time_steps_original);
    solve_rungekutta(n);
    double* y2 = copy_array(y_pointer, time_steps+1);

    set_time_steps(n1 * time_steps_original);
    solve_rungekutta(n);
    double* y1 = copy_array(y_pointer, time_steps+1);

    set_time_steps(time_steps_original);
    solve_rungekutta(n);
    double* y0 = copy_array(y_pointer, time_steps+1);

    for (int i=1; i<time_steps+1; i++){
        y_pointer[i] = y0[i] - pow(y1[n1*i] - y0[i], 2) / (y2[n2*i] - 2*y1[n1*i] + y0[i]);
    }
}

void IVP::solve_adams_exp_2(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }

    for (int i=(current_step-1); i<steps; i++){
        if (i<1){
            solve_rungekutta2(true, i+1);
        } else{
            t[i+1] = t[0] + (i+1)*h;
            y[i+1] = 3.*f_pointer(t[i], y[i]);
            y[i+1] += -f_pointer(t[i-1], y[i-1]);
            y[i+1] = y[i] + h/2. * y[i+1];
            f_evaluations += 1;
            f[i+1] = f_evaluations;
        }
    }
    reset_error_estimate();
}

void IVP::solve_adams_exp_3(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }

    for (int i=(current_step-1); i<steps; i++){
        if (i<2){
            solve_rungekutta3(true, i+1);
        } else{
            t[i+1] = t[0] + (i+1)*h;
            y[i+1] = 23.*f_pointer(t[i], y[i]);
            y[i+1] += -16.*f_pointer(t[i-1], y[i-1]);
            y[i+1] += 5.*f_pointer(t[i-2], y[i-2]);
            y[i+1] = y[i] + h/12. * y[i+1];
            f_evaluations += 1;
            f[i+1] = f_evaluations;
        }
    }
    reset_error_estimate();
}

void IVP::solve_adams_exp_4(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }

    for (int i=(current_step-1); i<steps; i++){
        if (i<3){
            solve_rungekutta4(true, i+1);
        } else{
            t[i+1] = t[0] + (i+1)*h;
            y[i+1] = 55.*f_pointer(t[i], y[i]);
            y[i+1] += -59.*f_pointer(t[i-1], y[i-1]);
            y[i+1] += 37.*f_pointer(t[i-2], y[i-2]);
            y[i+1] += -9.*f_pointer(t[i-3], y[i-3]);
            y[i+1] = y[i] + h/24. * y[i+1];
            f_evaluations += 1;
            f[i+1] = f_evaluations;
        }
    }
    reset_error_estimate();
}

void IVP::solve_adams_exp_5(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }

    for (int i=(current_step-1); i<steps; i++){
        if (i<4){
            solve_rungekutta5(true, i+1);
        } else{
            t[i+1] = t[0] + (i+1)*h;
            y[i+1] = 1901.*f_pointer(t[i], y[i]);
            y[i+1] += -2774.*f_pointer(t[i-1], y[i-1]);
            y[i+1] += 2616.*f_pointer(t[i-2], y[i-2]);
            y[i+1] += -1274.*f_pointer(t[i-3], y[i-3]);
            y[i+1] += 251.*f_pointer(t[i-4], y[i-4]);
            y[i+1] = y[i] + h/720. * y[i+1];
            f_evaluations += 1;
            f[i+1] = f_evaluations;
        }
    }
    reset_error_estimate();
}

void IVP::solve_adams_1(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }

    int count;
    double w_aprox, w_aprox_new, w_error;
    for (int i=(current_step-1); i<steps; i++){
        solve_adams_exp_2(true, i+1);

        w_aprox_new = y[i+1];
        w_error = 1;
        count = 0;
        while (w_error > conv_adams && count<max_iter_adams){
            w_aprox = w_aprox_new;

            w_aprox_new = f_pointer(t[i+1], w_aprox_new);
            w_aprox_new += f_pointer(t[i], y[i]);
            w_aprox_new = y[i] + h/2. * w_aprox_new;
            f_evaluations += 1;
            f[i+1] = f_evaluations;

            if (abs(w_aprox_new) > eps){
                w_error = abs((w_aprox - w_aprox_new)/w_aprox_new);
            } else{
                w_error = abs(w_aprox - w_aprox_new);
            }
            count++;
        }

        if (count > max_iter_adams && print_conv_adams){
            printf("Could not converge in time-step %d. Error: %g.\n",i+1,w_error);
        }

        y[i+1] = w_aprox_new;
    }
    reset_error_estimate();
}

void IVP::solve_adams_2(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }

    int count;
    double w_aprox, w_aprox_new, w_error;
    for (int i=(current_step-1); i<steps; i++){
        solve_adams_exp_3(true, i+1);
        if (i>0){
            w_aprox_new = y[i+1];
            w_error = 1;
            count = 0;
            while (w_error > conv_adams && count<max_iter_adams){
                w_aprox = w_aprox_new;

                w_aprox_new = 5.*f_pointer(t[i+1], w_aprox_new);
                w_aprox_new += 8.*f_pointer(t[i], y[i]);
                w_aprox_new += -f_pointer(t[i-1], y[i-1]);
                w_aprox_new = y[i] + h/12 * w_aprox_new;
                f_evaluations += 1;
                f[i+1] = f_evaluations;

                if (abs(w_aprox_new) > eps){
                    w_error = abs((w_aprox - w_aprox_new)/w_aprox_new);
                } else{
                    w_error = abs(w_aprox - w_aprox_new);
                }
                count++;
            }

            if (count > max_iter_adams && print_conv_adams){
                printf("Could not converge in time-step %d. Error: %g.\n",i+1,w_error);
            }

            y[i+1] = w_aprox_new;
        }
    }
    reset_error_estimate();
}

void IVP::solve_adams_3(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }

    int count;
    double w_aprox, w_aprox_new, w_error;
    for (int i=(current_step-1); i<steps; i++){
        solve_adams_exp_4(true, i+1);
        if (i>1){
            w_aprox_new = y[i+1];
            w_error = 1;
            count = 0;
            while (w_error > conv_adams && count<max_iter_adams){
                w_aprox = w_aprox_new;

                w_aprox_new = 9.*f_pointer(t[i+1], w_aprox_new);
                w_aprox_new += 19.*f_pointer(t[i], y[i]);
                w_aprox_new += -5.*f_pointer(t[i-1], y[i-1]);
                w_aprox_new += f_pointer(t[i-2], y[i-2]);
                w_aprox_new = y[i] + h/24. * w_aprox_new;
                f_evaluations += 1;
                f[i+1] = f_evaluations;

                if (abs(w_aprox_new) > eps){
                    w_error = abs((w_aprox - w_aprox_new)/w_aprox_new);
                } else{
                    w_error = abs(w_aprox - w_aprox_new);
                }
                count++;
            }

            if (count > max_iter_adams && print_conv_adams){
                printf("Could not converge in time-step %d. Error: %g.\n",i+1,w_error);
            }

            y[i+1] = w_aprox_new;
        }
    }
    reset_error_estimate();
}

void IVP::solve_adams_4(bool one_step, int current_step){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = h_current;
    double* t = t_pointer;
    double* y = y_pointer;
    int* f = evaluations_pointer;

    int steps = time_steps;
    if (one_step){
        steps = current_step;
    }

    int count;
    double w_aprox, w_aprox_new, w_error;
    for (int i=(current_step-1); i<steps; i++){
        solve_adams_exp_5(true, i+1);
        if (i>2){
            w_aprox_new = y[i+1];
            w_error = 1;
            count = 0;
            while (w_error > conv_adams && count<max_iter_adams){
                w_aprox = w_aprox_new;

                w_aprox_new = 251.*f_pointer(t[i+1], w_aprox_new);
                w_aprox_new += 646.*f_pointer(t[i], y[i]);
                w_aprox_new += -264.*f_pointer(t[i-1], y[i-1]);
                w_aprox_new += 106.*f_pointer(t[i-2], y[i-2]);
                w_aprox_new += -19.*f_pointer(t[i-3], y[i-3]);
                w_aprox_new = y[i] + h/720. * w_aprox_new;
                f_evaluations += 1;
                f[i+1] = f_evaluations;

                if (abs(w_aprox_new) > eps){
                    w_error = abs((w_aprox - w_aprox_new)/w_aprox_new);
                } else{
                    w_error = abs(w_aprox - w_aprox_new);
                }
                count++;
            }

            if (count > max_iter_adams && print_conv_adams){
                printf("Could not converge in time-step %d. Error: %g.\n",i+1,w_error);
            }

            y[i+1] = w_aprox_new;
        }
    }
    reset_error_estimate();
}

void  IVP::solve_adams(){
    initialize_problem();
    solve_adams(4, true);
}

void  IVP::solve_adams(int steps, bool implicit){
    initialize_problem();
    if (steps==1 && implicit){
        solve_adams_1(false,1);
    } else if (steps==2 && implicit){
        solve_adams_2(false,1);
    } else if (steps==2){
        solve_adams_exp_2(false,1);
    } else if (steps==3 && implicit){
        solve_adams_3(false,1);
    } else if (steps==3){
        solve_adams_exp_3(false,1);
    } else if (steps==4 && implicit){
        solve_adams_4(false,1);
    } else if (steps==4){
        solve_adams_exp_4(false,1);
    } else if (steps==5 && !implicit){
        solve_adams_exp_5(false,1);
    } else{
        printf("Invalid number of steps in Adams implementations!\n");
    }
}

void IVP::set_adams_convergence_limit(double value){
    conv_adams = value;
}

void IVP::set_adams_convergence_iterartions(int value){
    max_iter_adams = value;
}

void IVP::set_adams_convergence_print(bool value){
    print_conv_adams = value;
}

void IVP::solve_integral(){
    if (!has_solution()){
        printf("No solution found. Cannot continue.\n");
        return;
    }
    Integration integration;
    integration.add_points(t_pointer, y_pointer, time_steps+1);
    y_cum_pointer = integration.solve();
}

// ############# Tests #############

void tests_splines(){
    Spline spl;
    double* x = new double[3];
    double* y = new double[3];
    x[0] = 1;
    y[0] = 2;
    x[1] = 2;
    y[1] = 3;
    x[2] = 3;
    y[2] = 5;
    spl.add_points(x,y,3);
    printf("Added points\n");
    spl.build_natural();
    printf("Built natural spline\n");
    if (fabs(spl.get_y(2.4) - 3.704) > eps){
        printf("Error in y spline evaluation!\n");
    }
    if (fabs(spl.get_y_prime(2.4) - 1.98) > eps){
        printf("Error in y prime spline evaluation!\n");
    }
    if (fabs(spl.get_y_prime2(2.4) - 0.9) > eps){
        printf("Error in y prime2 spline evaluation!\n");
    }
    if (fabs(spl.get_int_y(2.4) - 3.7718999999999996) > eps){
        printf("Error in integral(y) spline evaluation!\n");
    }
    spl.add_y_prime(2,1);
    spl.build_fixed();
    printf("Built fixed spline\n");
    if (fabs(spl.get_y(2.4) - 3.824) > eps){
        printf("Error in y spline evaluation!\n");
    }
    if (fabs(spl.get_y_prime(2.4) - 2.38) > eps){
        printf("Error in y prime spline evaluation!\n");
    }
    if (fabs(spl.get_y_prime2(2.4) - 0.4) > eps){
        printf("Error in y prime2 spline evaluation!\n");
    }
    if (fabs(spl.get_int_y(2.4) - 3.8947333333333329) > eps){
        printf("Error in integral(y) spline evaluation!\n");
    }
}

double f_test_1(double t, double y){
    return y - t*t + 1;
}
double dfdt_test_1(double t, double y){
    return y - t*t + 1 - 2*t; // or delf_dely * f + delf_delt
}
double delfdely_test_1(double t, double y){
    return 1 + t*0 + y*0; // zero multiplication to avoid compilation warnings
}
double exact_test_1(double t){
    return pow(t+1, 2) - 0.5 * exp(t); // y(0)=0.5
}
double exact_int_test_1(double t){
    return t*(1. + t*(1. + t/3.) ) + 0.5 *(1. - exp(t));
}

double f_test_2(double t, double y){
    return -2.*y + 3.*exp(t);
}
double delfdely_test_2(double t, double y){
    return -2. + 0*t + 0*y;
}
double delfdelt_test_2(double t, double y){
    return 3.*exp(t) + 0*y;
}
double dfdt_test_2(double t, double y){
    return delfdely_test_2(t,y) * f_test_2(t,y) + delfdelt_test_2(t,y);
}
double exact_test_2(double t){
    return 2. * exp(-2.*t) + exp(t);  // y(0)=3
}
double exact_int_test_2(double t){
    return - exp(-2.*t) + exp(t);
}

double f_test_3(double t, double y){
    return 4*cos(t) - 8*sin(t) + 2*y;
}
double delfdely_test_3(double t, double y){
    return 2. + 0*t + 0*y;
}
double delfdelt_test_3(double t, double y){
    return -4*sin(t) - 8*cos(t) + 0*y;
}
double dfdt_test_3(double t, double y){
    return delfdely_test_3(t,y) * f_test_3(t,y) + delfdelt_test_3(t,y);
}
double exact_test_3(double t){
    return 4*sin(t) + 3*exp(2*t);  // y(0)=3
}
double exact_int_test_3(double t){
    return 4 - 4*cos(t) + 3*exp(t)*sinh(t);
}

double a = 3.;
double f_test_4(double t, double y){
    return -a*y + 0*t;
}
double delfdely_test_4(double t, double y){
    return -a + 0*t + 0*y;
}
double delfdelt_test_4(double t, double y){
    return 0. + 0*t + 0*y;
}
double dfdt_test_4(double t, double y){
    return delfdely_test_4(t,y) * f_test_4(t,y) + delfdelt_test_4(t,y);
}
double exact_test_4(double t){
    return 3.*exp(-a*t);  // y(0)=3
}
double exact_int_test_4(double t){
    return 3/a * (1- exp(-a*t));
}

void test_adams(IVP ivp, string problemname, string filename, FILE* evaluationsFile){
    const char* name = problemname.c_str();
    ivp.set_relative_error(false);

    ivp.set_adams_convergence_limit(1e-3);
    const int n = 10;

    printf(" Problem %s: Euler\n", name);
    ivp.set_time_steps(20*n);
    ivp.reset_f_evaluations();
    ivp.solve_euler();
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "Euler", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_euler.txt");

    printf(" Problem %s: Runge-Kutta 4\n", name);
    ivp.set_time_steps(5*n);
    ivp.reset_f_evaluations();
    ivp.solve_rungekutta(4);
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "RungeKutta4", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_rungekutta4.txt");

    printf(" Problem %s: Runge-Kutta 5\n", name);
    ivp.set_time_steps(n*5*2/3);
    ivp.reset_f_evaluations();
    ivp.solve_rungekutta(5);
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "RungeKutta5", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_rungekutta5.txt");

    printf(" Problem %s: Adams 2 step explicit\n", name);
    ivp.set_time_steps(20*n-1);
    ivp.reset_f_evaluations();
    ivp.solve_adams(2, false);
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "Adams2Exp", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_adams_2_exp.txt");

    printf(" Problem %s: Adams 3 step explicit\n", name);
    ivp.set_time_steps(20*n-4);
    ivp.reset_f_evaluations();
    ivp.solve_adams(3, false);
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "Adams3Exp", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_adams_3_exp.txt");

    printf(" Problem %s: Adams 4 step explicit\n", name);
    ivp.set_time_steps(20*n-9);
    ivp.reset_f_evaluations();
    ivp.solve_adams(4, false);
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "Adams4Exp", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_adams_4_exp.txt");

    printf(" Problem %s: Adams 5 step explicit\n", name);
    ivp.set_time_steps(20*n-20);
    ivp.reset_f_evaluations();
    ivp.solve_adams(5, false);
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "Adams5Exp", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_adams_5_exp.txt");

    printf(" Problem %s: Adams 1 step implicit\n", name);
    ivp.set_time_steps(10*n);
    ivp.reset_f_evaluations();
    ivp.solve_adams(1, true);
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "Adams1Imp", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_adams_1_imp.txt");

    printf(" Problem %s: Adams 2 step implicit\n", name);
    ivp.set_time_steps(10*n-1);
    ivp.reset_f_evaluations();
    ivp.solve_adams(2, true);
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "Adams2Imp", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_adams_2_imp.txt");

    printf(" Problem %s: Adams 3 step implicit\n", name);
    ivp.set_time_steps(10*n-3);
    ivp.reset_f_evaluations();
    ivp.solve_adams(3, true);
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "Adams3Imp", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_adams_3_imp.txt");

    printf(" Problem %s: Adams 4 step implicit\n", name);
    ivp.set_time_steps(10*n-8);
    ivp.reset_f_evaluations();
    ivp.solve_adams(4, true);
    printf("    'f' evaluations: %d\n", ivp.get_f_evaluations());
    fprintf(evaluationsFile, "%s\t%s\t%d\n", name, "Adams4Imp", ivp.get_f_evaluations());
    ivp.solve_integral();
    ivp.calculate_exact_error();
    // ivp.print_solution();
    ivp.print_solution(filename+"_adams_4_imp.txt");
}

void tests_adams(){

    FILE* outFile = fopen("evaluations.txt", "w");
    if (outFile == nullptr) {
        printf("Coud not create file: %s\n","evaluations.txt");
    }
    fprintf(outFile, "%s\t%s\t%s\n", "Problem", "Method", "Evaluations");

    printf("### Problem #1 ###\n");
    IVP test1;
    test1.set_f(f_test_1);
    test1.set_y_initial(0.5);
    test1.set_t_initial(0.);
    test1.set_t_end(2.);
    test1.set_exact(exact_test_1);
    test1.set_exact_cumulative(exact_int_test_1);
    test_adams(test1, "#1", "test1", outFile);

    printf("### Problem #2 ###\n");
    IVP test2;
    test2.set_f(f_test_2);
    test2.set_y_initial(3.);
    test2.set_t_initial(0.);
    test2.set_t_end(2.);
    test2.set_exact(exact_test_2);
    test2.set_exact_cumulative(exact_int_test_2);
    test_adams(test2, "#2", "test2", outFile);

    printf("### Problem #3 ###\n");
    IVP test3;
    test3.set_f(f_test_3);
    test3.set_y_initial(3.);
    test3.set_t_initial(0.);
    test3.set_t_end(2.);
    test3.set_exact(exact_test_3);
    test3.set_exact_cumulative(exact_int_test_3);
    test_adams(test3, "#3", "test3", outFile);

    printf("### Problem #4 ###\n");
    IVP test4;
    test4.set_f(f_test_4);
    test4.set_y_initial(3.);
    test4.set_t_initial(0.);
    test4.set_t_end(2.);
    test4.set_exact(exact_test_4);
    test4.set_exact_cumulative(exact_int_test_4);
    test_adams(test4, "#4", "test4", outFile);

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
double Fetkovich::get_f_pr(double We, double t) const{
    if (!has_f_pr()){
        printf("Reservoir pressure function not defined.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
    return f_pr_pointer(We,t);
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
            p_res[i] = get_f_pr(We[i], t[i]);
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


double f_pr_cte_230(double We, double t){
    return 230. + 0*We + 0*t;
}
double f_qw_pr_cte_230(double t){
    Fetkovich aqIVP1;
    aqIVP1.set_aquifer_initial_pore_volume(1*1E6);
    aqIVP1.set_aquifer_initial_pressure(250.);
    aqIVP1.set_aquifer_productivity_index(20.);
    aqIVP1.set_aquifer_total_compressibility(Newman_Consolidated_Sandstone(0.03));
    aqIVP1.set_reservoir_initial_pressure(230.);

    return aqIVP1.get_aquifer_flow_rate(t, 230.);
}

double f_pr_instant_res(double We, double t){
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

    double VoilIP_0 = Voil * Bob * (1 - Cob*(p_ini - Pb));
    double VwatIP_0 = Vwat * Bw;
    double Vpor_0 = VoilIP_0 + VwatIP_0;

    return (Voil * Bob * (1 + Cob*Pb) + (Vwat + We) * Bw - Vpor_0 * (1 - Cpor*(p_ini)) ) / (Vpor_0*Cpor + Voil * Bob * Cob) + 0*t;
}
double f_qw_instant_res(double t){
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

    double VoilIP_0 = Voil * Bob * (1 - Cob*(p_ini - Pb));
    double VwatIP_0 = Vwat * Bw;
    double Vpor_0 = VoilIP_0 + VwatIP_0;

    double J = 20.;
    double ct = Newman_Consolidated_Sandstone(0.03);
    double Wi = 1*1E6;
    double pi = 250.;
    double pr = 230.;

    return exp(-J * t *( 1. / (ct * Wi) + Bw / (Voil * Bob * Cob + Vpor_0 * Cpor))) * J * (pi - pr);
}
double f_qw_cumulative_res(double t){
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

    double VoilIP_0 = Voil * Bob * (1 - Cob*(p_ini - Pb));
    double VwatIP_0 = Vwat * Bw;
    double Vpor_0 = VoilIP_0 + VwatIP_0;

    double J = 20.;
    double ct = Newman_Consolidated_Sandstone(0.03);
    double Wi = 1*1E6;
    double pi = 250.;
    double pr = 230.;

    double alpha = ( 1. / (ct * Wi) + Bw / (Voil * Bob * Cob + Vpor_0 * Cpor));
    return (1 - exp(-J * t *alpha)) * (pi - pr) / alpha;
}
double f_instant_res(double t, double qw){
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

    double VoilIP_0 = Voil * Bob * (1 - Cob*(p_ini - Pb));
    double VwatIP_0 = Vwat * Bw;
    double Vpor_0 = VoilIP_0 + VwatIP_0;

    double J = 20.;
    double ct = Newman_Consolidated_Sandstone(0.03);
    double Wi = 1*1E6;
    // double pi = 250.;
    // double pr = 230.;

    return -J *( 1. / (ct * Wi) + Bw / (Voil * Bob * Cob + Vpor_0 * Cpor)) * qw + 0*t;
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
    aqFet.set_aquifer_initial_pressure(250.);
    aqFet.set_aquifer_productivity_index(20.);
    aqFet.set_aquifer_total_compressibility(Newman_Consolidated_Sandstone(0.03));
    aqFet.set_reservoir_initial_pressure(230.);
    aqFet.set_reservoir_pressure_function(f_pr_instant_res);
    aqFet.set_exact_water_flow_function(f_qw_instant_res);
    aqFet.set_exact_water_cumulative_function(f_qw_cumulative_res);

    aqFet.solve_aquifer_flow(200., 28*14/2);
    printf("    'f' evaluations: %d\n", aqFet.get_f_evaluations());
    aqFet.print_solution();
    aqFet.print_solution("aq1_fetkovich.txt");

    fprintf(outFile, "%s\t%s\t%d\n", "Aquifer#1", "Fetkovich", aqFet.get_f_evaluations());

    printf("Instant Pressure Equilibrium Reservoir as an IVP\n");
    IVP aqIVP1;
    aqIVP1.set_f(f_instant_res);
    aqIVP1.set_y_initial(400.);
    aqIVP1.set_t_initial(0.);
    aqIVP1.set_t_end(200.);
    aqIVP1.set_exact(f_qw_instant_res);
    aqIVP1.set_exact_cumulative(f_qw_cumulative_res);

    aqIVP1.set_adams_convergence_iterartions(200);
    test_adams(aqIVP1, "Aquifer#1", "aq1", outFile);
    fclose(outFile);

    int n_tests = 11;
    int* steps = new int[n_tests]{5, 10, 20, 50, 80, 100, 150, 200, 250, 500, 1000};

    double* error_list = new double;
    double error_end, error_max;
    int evaluations;

    printf("We Error Sensibility with Fetkovich\n");
    outFile = fopen("aq1_fetkovich_sens.txt", "w");
    printf("%10s\t%16s\t%16s\t%16s\n","Steps", "Evaluations", "ErrorEnd", "ErrorMax");
    fprintf(outFile,"%10s\t%16s\t%16s\t%16s\n","Steps", "Evaluations", "ErrorEnd", "ErrorMax");
    for (int i=0; i<10; i++){
        aqFet.reset_f_evaluations();
        aqFet.solve_aquifer_flow(200., steps[i]);
        evaluations = aqFet.get_f_evaluations();
        error_list = aqFet.get_result_cumulative_flow_error(true);
        error_end = error_list[steps[i]];
        error_max = max_array(error_list, steps[i]+1);
        printf("%10d\t%16d\t%16.10g\t%16.10g\n",steps[i], evaluations, error_end, error_max);
        fprintf(outFile,"%10d\t%16d\t%16.10g\t%16.10g\n",steps[i], evaluations, error_end, error_max);
    }
    fclose(outFile);

    string stringArray[9]={"rungekutta4","adams2exp","adams3exp","adams2impA","adams2impB","adams2impC","adams3impA","adams3impB","adams3impC"};

    aqIVP1.set_relative_error(true);
    aqIVP1.set_adams_convergence_iterartions(200);
    for (int j=0; j<9; j++){
        printf("We Error Sensibility with %s\n", stringArray[j].c_str());
        outFile = fopen(("aq1_" + stringArray[j] + "_sens.txt").c_str(), "w");
        printf("%10s\t%16s\t%16s\t%16s\n","Steps", "Evaluations", "ErrorEnd", "ErrorMax");
        fprintf(outFile,"%10s\t%16s\t%16s\t%16s\n","Steps", "Evaluations", "ErrorEnd", "ErrorMax");
        for (int i=0; i<10; i++){
            aqIVP1.set_time_steps(steps[i]);
            aqIVP1.reset_f_evaluations();
            switch (j) {
                case 0:
                    aqIVP1.solve_rungekutta(4);
                    break;
                case 1:
                    aqIVP1.solve_adams(2,false);
                    break;
                case 2:
                    aqIVP1.solve_adams(3,false);
                    break;
                case 3:
                    aqIVP1.set_adams_convergence_limit(1e-5);
                    aqIVP1.solve_adams(2,true);
                    break;
                case 4:
                    aqIVP1.set_adams_convergence_limit(1e-3);
                    aqIVP1.solve_adams(2,true);
                    break;
                case 5:
                    aqIVP1.set_adams_convergence_limit(1e-1);
                    aqIVP1.solve_adams(2,true);
                    break;
                case 6:
                    aqIVP1.set_adams_convergence_limit(1e-5);
                    aqIVP1.solve_adams(3,true);
                    break;
                case 7:
                    aqIVP1.set_adams_convergence_limit(1e-3);
                    aqIVP1.solve_adams(3,true);
                    break;
                case 8:
                    aqIVP1.set_adams_convergence_limit(1e-1);
                    aqIVP1.solve_adams(3,true);
                    break;
                default:
                    printf("Undifined!\n");
                    return;
            }
            evaluations = aqIVP1.get_f_evaluations();
            aqIVP1.solve_integral();
            aqIVP1.calculate_exact_error();
            error_list = aqIVP1.get_y_cumulative_error();
            error_end = error_list[steps[i]];
            error_max = max_array(error_list, steps[i]+1);
            printf("%10d\t%16d\t%16.10g\t%16.10g\n",steps[i], evaluations, error_end, error_max);
            fprintf(outFile,"%10d\t%16d\t%16.10g\t%16.10g\n",steps[i], evaluations, error_end, error_max);
        }
        fclose(outFile);
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
    // tests_splines();
    // tests_integration();
    tests_adams();
    // Fetkovich_tests();
}