#include <iostream>
#include <cmath> 
#include <vector>
using namespace std;

const double eps = 1E-12;

class IVP{
    public:
        IVP();
        void set_f(double (*f)(double, double));
        void set_exact(double (*f)(double));
        void set_y_initial(double y);
        void set_t_initial(double t);
        void set_t_end(double t);
        void set_time_steps(int n);

        double get_f(double t, double y);
        double get_exact(double t);

        bool has_f();
        bool has_exact();
        bool has_solution();

        void solve_euler();

        void print_solution();
    private:
        double (*f_pointer)(double, double);
        double (*exact_pointer)(double);
        double y_initial;
        double t_initial;
        double t_end;
        int time_steps;
        double *t_pointer;
        double *y_pointer;
        bool relative_error;
};

IVP::IVP(): f_pointer(nullptr), exact_pointer(nullptr), y_initial(0), t_initial(0), t_end(1), time_steps(100), t_pointer(nullptr), y_pointer(nullptr), relative_error(false) {
};
void IVP::set_f(double (*f)(double, double)){
    f_pointer = f;
}
void IVP::set_exact(double (*f)(double)){
    exact_pointer = f;
}
void IVP::set_y_initial(double y){
    y_initial = y;
}
void IVP::set_t_initial(double t){
    t_initial = t;
    if (abs(t-t_end)<eps){
        t_end = t_initial + 1;
        printf("Cannot set initial time equal to end time. Setting end time to initial time + 1: %g\n", t_end);
    }
}
void IVP::set_t_end(double t){
    if (abs(t-t_initial)<eps){
        printf("Cannot set end time equal to initial time. Keeping previous value: %g\n", t_end);
    } else{
        t_end = t;
    }
}
void IVP::set_time_steps(int n){
    if (n<1){
        printf("Cannot set less than one time-step. Keeping previous value: %i\n", time_steps);
    } else{
        time_steps = n;
    }
}

double IVP::get_f(double t, double y){
    return f_pointer(t,y);
}
double IVP::get_exact(double t){
    return exact_pointer(t);
}

bool IVP::has_f(){
    return f_pointer != nullptr;
}
bool IVP::has_exact(){
    return exact_pointer != nullptr;
}
bool IVP::has_solution(){
    return t_pointer != nullptr && y_pointer != nullptr;
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
}

void IVP::print_solution(){
    if (!has_solution()){
        printf("There is no solution to print.\n");
        return;
    }
    if (!has_exact()){
        printf("%5s\t%16s\t%16s\n","#","t","y");
        for (int i=0; i<time_steps+1; i++){
            printf("%5i\t%16.10g\t%16.10g\n", i, t_pointer[i], y_pointer[i]);
        }
    } else{
        double y_calc;
        double y_error;

        printf("%5s\t%16s\t%16s\t%16s\t%16s\n","#","t","y aprox.","y real","y error");
        for (int i=0; i<time_steps+1; i++){
            y_calc = get_exact(t_pointer[i]);
            if (relative_error && abs(y_calc) > eps){
                y_error = abs((y_calc - y_pointer[i])/y_calc);
            } else{
                y_error = abs(y_calc - y_pointer[i]);
            }
            printf("%5i\t%16.10g\t%16.10g\t%16.10g\t%16.10g\n", i, t_pointer[i], y_pointer[i],y_calc,y_error);
        }
    }
}

double f_test_1(double t, double y){
    return y - t*t + 1;
}
double exact_test_1(double t){
    return (t+1)*(t+1) - 0.5 * exp(t);
}

void tests_euler(){
    printf("### Verify 'error' catching ###\n");
    IVP test_errors;
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
    test1.print_solution();
    printf("\n");



    if (test1.has_f()){
        printf("Tem f\n");
    } else{
        printf("Nao tem f\n");
    }
}

int main(){
    tests_euler();
}