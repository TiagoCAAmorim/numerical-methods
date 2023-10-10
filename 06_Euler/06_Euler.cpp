#include <iostream>
#include <cmath> 
using namespace std;

const double eps = 1E-12;

class IVP{
    public:
        IVP();
        
        void set_f(double (*f)(double, double)){
            f_pointer = f;
        }
        void set_solution(double (*f)(double)){
            solution_pointer = f;
        }
        void set_y_initial(double y){
            y_initial = y;
        }
        void set_t_initial(double t){
            t_initial = t;
        }
        void set_t_end(double t){
            if (abs(t-t_initial)<eps){
                printf("Cannot set end time equal to initial time. Try again.\n");
            } else{
                t_end = t;
            }
        }
        void set_time_steps(int n){
            if (n<1){
                printf("Cannot set less than one time-step.\n");
            } else{
                time_steps = n;
            }
        }

        get_f(double t, double y){
            return f_pointer(t,y);
        }
        get_solution(double t){
            return solution_pointer(t);
        }

        bool has_f(){
            return f_pointer != nullptr;
        }
        bool has_solution(){
            return solution_pointer != nullptr;
        }
    private:
        double (*f_pointer)(double, double);
        double (*solution_pointer)(double);
        double y_initial;
        double t_initial;
        double t_end;
        int time_steps;

};

IVP::IVP(): f_pointer(nullptr), solution_pointer(nullptr), y_initial(0), t_initial(0), t_end(1), time_steps(100) {
};


void euler(double (*f)(double, double), double t_initial, double t_end, double y_initial, int time_steps, double *t,double *y){
    double h = (t_end - t_initial) / time_steps;
    t[0] = t_initial;
    y[0] = y_initial;
    
    for (int i=1; i<time_steps+1; i++){
        y[i] = y[i-1] + h * f(t[i-1], y[i-1]);
        t[i] = t[0] + i*h;
    }
}

void print_vectors(int n, double *t, double *y){
    printf("%5s\t%16s\t%16s\n","#","t","y");
    for (int i=0; i<n; i++){
        printf("%5i\t%16.10g\t%16.10g\n", i, t[i], y[i]);
    }
}

void print_vectors(int n, double *t, double *y, double (*y_real)(double), bool relative_error=false){
    double y_calc;
    double y_error;

    printf("%5s\t%16s\t%16s\t%16s\t%16s\n","#","t","y aprox.","y real","y error");
    for (int i=0; i<n; i++){
        y_calc = y_real(t[i]);
        if (relative_error && abs(y_calc) > eps){
            y_error = abs((y_calc - y[i])/y_calc);
        } else{
            y_error = abs(y_calc - y[i]);
        }
        printf("%5i\t%16.10g\t%16.10g\t%16.10g\t%16.10g\n", i, t[i], y[i],y_calc,y_error);
    }
}




double f_test_1(double t, double y){
    return y - t*t + 1;
}

void tests_euler(){
    IVP test1;
    if (test1.has_f()){
        printf("Tem f\n");
    } else{
        printf("Nao tem f\n");
    }
}

int main(){
    tests_euler();
}