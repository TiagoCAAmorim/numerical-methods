#include <iostream>
#include <cmath>
#include <limits>
#include <cstdio>
#include <string>
using namespace std;

const double eps = 1E-12;

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

        // void copy_array(double *original, double *new_vector, int number_points) const;
        double* copy_array(double *original, int number_points) const;

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

double* Spline::copy_array(double *original, int number_points) const{
    double* out = new double[number_points];
    for(int i=0; i<number_points; i++){
        out[i] = original[i];
    }
    return out;
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

    if (fabs(b[0]  - y_prime_pointer[0]) > eps){
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


class IVP{
    public:
        IVP();
        void set_f(double (*f)(double, double));
        void set_exact(double (*f)(double));
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

        double get_f(double t, double y) const;
        double get_dfdt(double t, double y) const;
        double get_delfdely(double t, double y) const;
        double get_exact(double t) const;
        double get_Lipschitz() const;
        double get_max_dfdt() const;

        void calculate_exact_error();
        void print_solution();
        void print_solution(string filename);

        void solve_euler();
        void estimate_error_euler();
    private:
        void reset_results();
        void reset_error_estimate();

        bool has_f() const;
        bool has_exact() const;
        bool has_solution() const;
        bool has_exact_error() const;
        bool has_estimated_error() const;
        bool has_dfdt() const;
        bool has_delfdely() const;
        bool has_Lipschitz() const;
        bool has_max_dfdt() const;

        double (*f_pointer)(double, double);
        double (*exact_pointer)(double);
        double (*dfdt_pointer)(double, double);
        double (*delfdely_pointer)(double, double);
        double y_initial;
        double t_initial;
        double t_end;
        int time_steps;
        double *t_pointer;
        double *y_pointer;
        double *y_exact_pointer;
        double *y_error_pointer;
        double *y_error_limit_pointer;
        int extra_search_points;
        bool relative_error;
        double Lipschitz;
        bool calculated_L;
        double max_dfdt;
        bool calculated_M;

};

IVP::IVP():
    f_pointer(nullptr), exact_pointer(nullptr),
    dfdt_pointer(nullptr),
    delfdely_pointer(nullptr),
    y_initial(0), t_initial(0),
    t_end(1), time_steps(100),
    t_pointer(nullptr), y_pointer(nullptr),
    y_exact_pointer(nullptr), y_error_pointer(nullptr),
    y_error_limit_pointer(nullptr),
    extra_search_points(20), relative_error(false),
    Lipschitz(0), calculated_L(false),
    max_dfdt(0), calculated_M(false)
    {};

void IVP::reset_results(){
    t_pointer = nullptr;
    y_pointer = nullptr;
    calculated_L = false;
    calculated_M = false;
    reset_error_estimate();
}
void IVP::reset_error_estimate(){
    y_exact_pointer = nullptr;
    y_error_pointer = nullptr;
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
    relative_error = relative;
    y_exact_pointer = nullptr;
    y_error_pointer = nullptr;
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

double IVP::get_f(double t, double y) const{
    if (!has_f()){
        printf("Function f not defined.\n");
        return std::numeric_limits<double>::quiet_NaN();
    }
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
bool IVP::has_solution() const{
    return t_pointer != nullptr && y_pointer != nullptr;
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
}

void IVP::print_solution(){
    if (!has_solution()){
        printf("There is no solution to print.\n");
        return;
    }

    printf("%5s\t%16s\t%16s","#","t","y_aprox");
    if (has_estimated_error()){
        printf("\t%16s","error_limit");
    }
    if (has_exact_error()){
        printf("\t%16s\t%16s","y_exact","exact_error");
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

    fprintf(outFile,"%5s\t%16s\t%16s","#","t","y_aprox");
    if (has_estimated_error()){
        fprintf(outFile,"\t%16s","error_limit");
    }
    if (has_exact_error()){
        fprintf(outFile,"\t%16s\t%16s","y_exact","exact_error");
    }
    fprintf(outFile,"\n");

    for (int i=0; i<time_steps+1; i++){
        fprintf(outFile,"%5i\t%16.10g\t%16.10g", i, t_pointer[i], y_pointer[i]);
        if (has_estimated_error()){
            fprintf(outFile,"\t%16.10g", y_error_limit_pointer[i]);
        }
        if (has_exact_error()){
            fprintf(outFile,"\t%16.10g\t%16.10g", y_exact_pointer[i],y_error_pointer[i]);
        }
        fprintf(outFile,"\n");
    }
    fclose(outFile);
}

void IVP::solve_euler(){
    if (!has_f()){
        printf("Function f(t,y) not defined. Cannot continue.\n");
        return;
    }
    double h = (t_end - t_initial) / time_steps;
    double* t = new double[time_steps+1];
    double* y = new double[time_steps+1];

    t[0] = t_initial;
    y[0] = y_initial;

    for (int i=1; i<time_steps+1; i++){
        y[i] = y[i-1] + h * f_pointer(t[i-1], y[i-1]);
        t[i] = t[0] + i*h;
    }

    t_pointer = t;
    y_pointer = y;
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
    return y - t*t + 1 - 2*t;
}
double delfdely_test_1(double t, double y){
    return 1 + t*0 + y*0; // zero multiplication to avoid compilation warnings
}
double exact_test_1(double t){
    return (t+1)*(t+1) - 0.5 * exp(t);
}

void tests_euler(){
    printf("### Verify error catching ###\n");
    IVP test_errors;
    test_errors.get_f(0,1);
    test_errors.get_exact(0);
    test_errors.get_dfdt(0,1);
    test_errors.get_delfdely(0,1);
    test_errors.get_Lipschitz();
    test_errors.estimate_Lipschitz();
    test_errors.solve_euler();
    test_errors.set_t_end(0.);
    test_errors.set_t_initial(1.);
    test_errors.set_time_steps(0);
    test_errors.print_solution();
    printf("\n");

    printf("### Problem #1 ###\n");
    IVP test1;
    test1.set_f(f_test_1);
    test1.set_y_initial(0.5);
    test1.set_t_initial(0.);
    test1.set_t_end(2.);
    test1.set_time_steps(10);
    test1.set_exact(exact_test_1);
    test1.solve_euler();
    test1.set_relative_error(false);
    test1.calculate_exact_error();

    printf("  With aproximate L and M\n");
    test1.set_extra_search_points(3);
    test1.set_delfdely(delfdely_test_1);
    test1.estimate_Lipschitz();
    test1.set_dfdt(dfdt_test_1);
    test1.estimate_max_dfdt();
    printf("L = %g\n",test1.get_Lipschitz());
    printf("M = %g\n",test1.get_max_dfdt());
    test1.estimate_error_euler();
    test1.print_solution();
    test1.print_solution("test1.txt");

    printf("  With exact L and M\n");
    test1.set_Lipschitz_L(1);
    test1.set_max_dfdt(0.5*exp(2)-2);
    printf("L = %g\n",test1.get_Lipschitz());
    printf("M = %g\n",test1.get_max_dfdt());
    test1.estimate_error_euler();
    test1.print_solution();

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

        double get_aquifer_flow_rate(double t, double pr) const;  // m3/d
        double get_aquifer_delta_cumulative_flow(double dt, double paq_avg, double pr_avg) const;  // m3

        void solve_aquifer_flow(double t_end, int steps);
        double* get_result_times() const; // d
        double* get_result_water_flow() const; // m3/d
        double* get_result_cumulative_flow_rate() const; // m3
        double* get_result_aquifer_pressure() const; // bar
        double* get_result_reservoir_pressure() const; // bar

        void print_solution(string filename);
    private:
        double get_aquifer_pressure(double We) const;   // bar

        bool has_f_pr() const;
        double get_f_pr(double We, double t) const;

        bool has_exact() const;
        double get_exact(double t) const;

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

        // double Qw;
        // double We;
        // double Pa;
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
    f_exact_pointer(nullptr)
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

double Fetkovich::get_aquifer_flow_rate(double t, double pr) const{
    return exp(-J * t / (ct * Wi)) * J * (pi - pr);
}
double Fetkovich::get_aquifer_delta_cumulative_flow(double dt, double paq_avg, double pr_avg) const{
    return (1 - exp(-J * dt / (ct * Wi))) * ct * Wi * (paq_avg - pr_avg);
}
double Fetkovich::get_aquifer_pressure(double We) const{
    return pi - We / (ct * Wi);
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
        while (abs(pr_avg_old - pr_avg) > eps && count<=20){
            pr_avg_old = pr_avg;
            dWe = get_aquifer_delta_cumulative_flow(dt, p_aq[i-1], pr_avg);
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
double* Fetkovich::get_result_cumulative_flow_rate() const{
    return We_pointer;
}
double* Fetkovich::get_result_aquifer_pressure() const{
    return p_aquifer_pointer;
}
double* Fetkovich::get_result_reservoir_pressure() const{
    return p_reservoir_pointer;
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

    fprintf(outFile,"%5s\t%16s","#","Time");
    fprintf(outFile,"\t%16s","Res.Pres.");
    fprintf(outFile,"\t%16s","Aq.Pres.");
    fprintf(outFile,"\t%16s","Wat.Flow");
    fprintf(outFile,"\t%16s","Cumulative");
    if (has_exact()){
        fprintf(outFile,"\t%16s","ExactQw");
        fprintf(outFile,"\t%16s","Error");
        fprintf(outFile,"\t%16s","ErrorRel");
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
        fprintf(outFile,"\n");
    }
    fclose(outFile);
}



double Newmam_Consolidated_Sandstone(double por){
    if (por >= 1.){
        por = por / 100.;
    }
    return 97.32E-6 / (1 + 55.8721 * pow(por, 1.428586) ) * 14.50377; // 1/bar
}
double Newmam_Limestone(double por){
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
    Fetkovich aq1;
    aq1.set_aquifer_initial_pore_volume(1*1E6);
    aq1.set_aquifer_initial_pressure(250.);
    aq1.set_aquifer_productivity_index(20.);
    aq1.set_aquifer_total_compressibility(Newmam_Consolidated_Sandstone(0.03));
    aq1.set_reservoir_initial_pressure(230.);

    return aq1.get_aquifer_flow_rate(t, 230.);
}

double f_pr_instant_res(double We, double t){
    double api = 25.;
    double dg = 0.6;
    double rgo = 60.;
    double temp = 65.;
    double Bob = Standing_bo_bubble(api, dg, rgo, temp);
    double Cob = Standing_co_bubble(api, dg, rgo, temp);
    double Pb  = Standing_p_bubble(api, dg, rgo, temp);
    double Cpor = Newmam_Consolidated_Sandstone(0.2);
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
    double Cpor = Newmam_Consolidated_Sandstone(0.2);
    double Bw = 1.;
    double Voil = 0.8 * 1*1E6;
    double Vwat = 0.2 * 1*1E6;
    double p_ini = 230.;

    double VoilIP_0 = Voil * Bob * (1 - Cob*(p_ini - Pb));
    double VwatIP_0 = Vwat * Bw;
    double Vpor_0 = VoilIP_0 + VwatIP_0;

    double J = 20.;
    double ct = Newmam_Consolidated_Sandstone(0.03);
    double Wi = 1*1E6;
    double pi = 250.;
    double pr = 230.;

    return exp(-J * t *( 1. / (ct * Wi) + Bw / (Voil * Bob * Cob + Vpor_0 * Cpor))) * J * (pi - pr);
}

void Fetkovich_tests(){
    Fetkovich aq1;
    aq1.set_aquifer_initial_pore_volume(1*1E6);
    aq1.set_aquifer_initial_pressure(250.);
    aq1.set_aquifer_productivity_index(20.);
    aq1.set_aquifer_total_compressibility(Newmam_Consolidated_Sandstone(0.03));
    aq1.set_reservoir_initial_pressure(230.);
    aq1.set_reservoir_pressure_function(f_pr_cte_230);
    aq1.set_exact_water_flow_function(f_qw_pr_cte_230);

    aq1.solve_aquifer_flow(300., 10);
    aq1.print_solution("aq1.txt");

    Fetkovich aq2;
    aq2.set_aquifer_initial_pore_volume(1*1E6);
    aq2.set_aquifer_initial_pressure(250.);
    aq2.set_aquifer_productivity_index(20.);
    aq2.set_aquifer_total_compressibility(Newmam_Consolidated_Sandstone(0.03));
    aq2.set_reservoir_initial_pressure(230.);
    aq2.set_reservoir_pressure_function(f_pr_instant_res);
    aq2.set_exact_water_flow_function(f_qw_instant_res);

    aq2.solve_aquifer_flow(300., 10);
    aq2.print_solution("aq2.txt");
}





int main(){
    // tests_splines();
    // tests_euler();
    Fetkovich_tests();
}