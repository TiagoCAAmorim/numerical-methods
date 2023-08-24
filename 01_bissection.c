/*
    Implementation of the bissection method to find root of 1D functions    
*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

const double epsilon = 1E-7;
const double pi = 3.14159265358979323846;
// const double error_number = -9999.8999;

bool is_root(double fx, double x, int i, bool debug){
    if( fabs(fx) < epsilon){
        if( debug)
            printf("'Root' found after %i interations.\nFound f(x) = %g at x = %g.\n", i, fx, x);
        return true;
    }
    return false;
}

double calculate_convergence(double x_current, double x_previous, bool relative_convergence){
    double reference = 1;
    if( relative_convergence)
        if( fabs(x_current) > epsilon){
            reference = x_current;
        } else if( fabs(x_previous) > epsilon) {
            reference = x_previous;
        } else{
            return 0.;
        }
    return fabs( (x_current - x_previous) / reference );
}

double find_root_bissection_debug(double (*func)(double), double x_a, double x_b, double convergence_tol, bool relative_convergence, int max_interations, bool debug) {
    double fx_a, fx_b, fx_mean;
    double x_mean, x_mean_previous;
    double convergence;

    if( debug)
        printf("### Bissection Method ###\n");

    fx_a = func(x_a);
    if( is_root(fx_a, x_a, 0, debug))
        return x_a;

    fx_b = func(x_b);
    if( is_root(fx_b, x_b, 0, debug))
        return x_b;

    if( signbit(fx_a) == signbit(fx_b) ){
        printf("Function has same sign in limits:\n f(%g)=%g\tf(%g)=%g\nCannot proceed.\n",x_a,fx_a,x_b,fx_b);
        assert(false);
    }

    x_mean_previous = x_a;
    for(int i=1 ; i<=max_interations ; i++){
        x_mean = x_a + (x_b - x_a)/2;
        convergence = calculate_convergence(x_mean, x_mean_previous, relative_convergence);
        x_mean_previous = x_mean;
        
        if( convergence < convergence_tol){
            if( debug)
                printf("Reached convergence after %i interations at x = %g.\n", i, x_mean);
            return x_mean;
        }

        fx_mean = func(x_mean);
        if( debug)
            printf("%3i.  f(%-10g) = %-10g\tf(%-10g) = %-10g\tConvergence=%-10g\n",i,x_a,fx_a,x_b,fx_b, convergence);
        
        if( is_root(fx_mean, x_mean, i, debug))
            return x_mean;

        if( signbit(fx_a) == signbit(fx_mean)){
            x_a = x_mean;
            fx_a = fx_mean;
        } else{
            x_b = x_mean;
            fx_b = fx_mean;
        }
    }
    if( debug)
        printf("Convergence was not reached after %i interations.\nFound f(x) = %g at x = %g.\n", max_interations, fx_mean, x_mean);
    return x_mean;
}

double find_root_bissection(double (*func)(double), double x_a, double x_b, double convergence_tol, bool relative_convergence, int max_interations) {
    return find_root_bissection_debug(func, x_a, x_b, convergence_tol, relative_convergence, max_interations, false);
}

double min(double a, double b){
    return a<b? a: b;
}

double max(double a, double b){
    return a>b? a: b;
}

int estimate_bissection_interations(double x_a, double x_b, double x_root, double convergence_tol, bool relative_convergence){
    double x_reference = 1;
    double n;

    if( relative_convergence)
        if(x_root<epsilon || x_root < min(x_a,x_b) || x_root > max(x_a,x_b)){
            x_reference = min(fabs(x_a), fabs(x_b));
            if( x_reference < epsilon)
                x_reference = max(fabs(x_a), fabs(x_b));
            if( x_reference < epsilon)
                return 0;
        } else{
            x_reference = x_root;
        }
    n = log(fabs(x_b - x_a) / fabs(x_reference) / convergence_tol)/log(2);
    return round(n+0.5);
}

void test_bissection(double (*func)(double), double x_root, double x_a, double x_b, double convergence_tol, bool relative_convergence, int max_interations, bool debug, char *message){
    printf("Test Bissection Method: %s\n", message);
    int interations = estimate_bissection_interations(x_a, x_b, x_root, convergence_tol, false);
    printf("  Estimated number of interations: %i\n", interations);
    double x_root_bissection = find_root_bissection_debug(func, x_a, x_b, convergence_tol, relative_convergence, max_interations, debug); 
    double convergence = calculate_convergence(x_root, x_root_bissection, false);
    printf("  Root=%g. Found x=%g.\nError on root: %g.\n\n", x_root, x_root_bissection, convergence);
}

// Root at x=0.3
double f_linear(double x){
    return -x*3 + 0.9;
}

// Roots at x=-0.5 and -0.1
double f_quadratic(double x){
    return 5*x*x + 3*x + 0.25;
}

// Root at x= pi/2 (+ n*2pi)
double f_cos(double x){
    return cos(x);
}

// Root at x= 3/4*pi (+ n*2pi)
double f_trigonometric(double x){
    return cos(x) + sin(x);
}

int main(){

    test_bissection(f_linear, 0.3, 0, 1, 0.001, false, 50, true, "Linear function: ");
    test_bissection(f_linear, 0.3, 0.3, 1, 0.001, false, 50, false, "Linear function with root in x_a: ");
    test_bissection(f_linear, 0.3, 0, 0.3, 0.001, false, 50, false, "Linear function with root in x_b: ");
    // test_bissection(f_linear, 0.3, 1, 2.3, 0.001, false, 50, false, "Linear function with error in [x_a,x_b]: ");
    // test_bissection(f_linear, 0.3, -2, 0., 0.001, false, 50, false, "Linear function with error in [x_a,x_b]: ");
    
    test_bissection(f_quadratic, -0.1, -0.25, 1, 0.001, false, 50, true, "Quadratic function, root#1: ");
    test_bissection(f_quadratic, -0.5, -0.25, -1, 0.001, false, 50, true, "Quadratic function, root#2: ");
    
    test_bissection(f_cos, pi/2, 0, 2, 0.001, false, 50, true, "Cossine function: ");
    test_bissection(f_trigonometric, 3./4*pi, 0, 5, 0.001, false, 50, true, "Trigonometric function: ");

    return 0;
}