/*
    Implementation of the bissection method to find root of 1D functions    
*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

const float epsilon = 1E-7;

bool is_root(float fx, float x, int i, bool debug){
    if( fabs(fx) < epsilon){
        if( debug)
            printf("'Root' found after %i interations.\nFound f(x) = %g at x = %g.\n", i, fx, x);
        return true;
    }
    return false;
}

float calculate_convergence(float x_current, float x_previous){
    if( fabs(x_current) > epsilon){
        return fabs( (x_current - x_previous) / x_current );
    } else if( fabs(x_previous) > epsilon) {
        return fabs( (x_current - x_previous) / x_previous );
    } else{
        return 0.;
    }
}

float find_root_bissection_debug(float (*func)(float), float x_a, float x_b, float convergence_tol, int max_interations, bool debug) {
    float fx_a, fx_b, fx_mean;
    float x_mean, x_mean_previous;
    float convergence;

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
        convergence = calculate_convergence(x_mean, x_mean_previous);
        x_mean_previous = x_mean;
        
        if( convergence < convergence_tol){
            if( debug)
                printf("Reached convergence after %i interations at x = %g.\n", i, x_mean);
            return x_mean;
        }

        fx_mean = func(x_mean);
        if( debug)
            printf("%i.\tf(%g)=%g\tf(%g)=%g\n",i,x_a,fx_a,x_b,fx_b);
        
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

float find_root_bissection(float (*func)(float), float x_a, float x_b, float convergence_tol, int max_interations) {
    return find_root_bissection_debug(func, x_a, x_b, convergence_tol, max_interations, false);
}

// Root at x=0.3
float f_linear(float x){
    return -x*3 + 0.9;
}

// Roots at x=-0.5 and -0.1
float f_quadratic(float x){
    return 5*x*x + 3*x + 0.25;
}

void test_bissection(float (*func)(float), float x_root, float x_a, float x_b, float convergence_tol, int max_interations, bool debug, char *message){
    float x_root_bissection = find_root_bissection_debug(func, x_a, x_b, convergence_tol, max_interations, debug); 
    float convergence = calculate_convergence(x_root, x_root_bissection);
    printf("%sRoot=%g\tFound x=%g\nError on root: %g\n\n", message, x_root, x_root_bissection, convergence);
}

int main(){

    test_bissection(f_linear, 0.3, 0, 1, 0.001, 50, true, "Linear function: ");
    test_bissection(f_linear, 0.3, 0.3, 1, 0.001, 50, false, "Linear function with root in x_a: ");
    test_bissection(f_linear, 0.3, 0, 0.3, 0.001, 50, false, "Linear function with root in x_b: ");
    // test_bissection(f_linear, 0.3, 1, 2.3, 0.001, 50, false, "Linear function with error in [x_a,x_b]: ");
    // test_bissection(f_linear, 0.3, -2, 0., 0.001, 50, false, "Linear function with error in [x_a,x_b]: ");
    
    test_bissection(f_quadratic, -0.1, -0.25, 1, 0.001, 50, true, "Quadratic function, root#1: ");
    test_bissection(f_quadratic, -0.5, -0.25, -1, 0.001, 50, true, "Quadratic function, root#2: ");

    return 0;
}