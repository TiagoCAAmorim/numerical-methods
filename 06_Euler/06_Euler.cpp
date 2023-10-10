#include <iostream>
#include <cmath> 
#include <limits>
using namespace std;

const double eps = 1E-12;

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

        void solve_euler();
        void estimate_error_euler();
    private:
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
};

IVP::IVP(): f_pointer(nullptr), exact_pointer(nullptr),
            dfdt_pointer(nullptr), 
            delfdely_pointer(nullptr), 
            y_initial(0), t_initial(0), 
            t_end(1), time_steps(100), 
            t_pointer(nullptr), y_pointer(nullptr), 
            y_exact_pointer(nullptr), y_error_pointer(nullptr), 
            extra_search_points(20), relative_error(false),
            y_error_limit_pointer(nullptr),
            Lipschitz(0), calculated_L(false),
            max_dfdt(0), calculated_M(false) {
};
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

    if (extra_search_points > 0){
        //build spline
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
            // y = spl.get_value(t);
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

    if (extra_search_points > 0){
        //build spline
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
            // y = spl.get_value(t);
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
        if (relative_error){
            if (abs(y_exact[i]) > eps){
                y_error[i] = abs((y_exact[i] - y_pointer[i])/y_exact[i]);
            } else{
                y_error[i] = std::numeric_limits<double>::quiet_NaN();
            }
        } else{
            y_error[i] = abs(y_exact[i] - y_pointer[i]);
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
    }

    y_error_limit_pointer = error_limit;
}




double f_test_1(double t, double y){
    return y - t*t + 1;
}
double dfdt_test_1(double t, double y){
    return y - t*t + 1 - 2*t;
}
double delfdely_test_1(double t, double y){
    return 1;
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
    test1.calculate_exact_error();

    printf("  With aproximate L and M\n");
    test1.set_delfdely(delfdely_test_1);
    test1.estimate_Lipschitz();
    test1.set_dfdt(dfdt_test_1);
    test1.estimate_max_dfdt();
    printf("L = %g\n",test1.get_Lipschitz());
    printf("M = %g\n",test1.get_max_dfdt());
    test1.estimate_error_euler();
    test1.print_solution();
    printf("  With exact L and M\n");
    test1.set_Lipschitz_L(1);
    test1.set_max_dfdt(0.5*exp(2)-2);
    printf("L = %g\n",test1.get_Lipschitz());
    printf("M = %g\n",test1.get_max_dfdt());
    test1.estimate_error_euler();
    test1.print_solution();

}

int main(){
    tests_euler();
}