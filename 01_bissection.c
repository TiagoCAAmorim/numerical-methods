/*
    Implementation of the bissection method to find root of 1D functions    
*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

const double epsilon = 1E-7;
const double pi = 3.14159265358979323846;

const int default_interations = 100;
const double default_convergence = 1E-4;

double min(double a, double b){
    return a<b? a: b;
}

double max(double a, double b){
    return a>b? a: b;
}

bool is_root(double fx, double x, int i, bool debug){
    if( fabs(fx) < epsilon){
        if( debug){
            printf("Early exit after %i interations: f(%g) = %g.\n", i, x, fx);
        }
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

void print_error(char *message, double true_value, double calculated_value){
    double relative_error;
    relative_error = calculate_convergence(true_value, calculated_value, true);
    printf("%s: True=%g Calculated=%g Error=%g\n", message, true_value, calculated_value, relative_error);
}

double find_root_bissection_debug(double (*func)(double), double x_a, double x_b, double convergence_tol, bool relative_convergence, int max_interations, bool debug) {
    double fx_a, fx_b, fx_mean;
    double x_mean, x_mean_previous;
    double convergence;

    fx_a = func(x_a);
    if( is_root(fx_a, x_a, 0, debug))
        return x_a;

    fx_b = func(x_b);
    if( is_root(fx_b, x_b, 0, debug))
        return x_b;

    if( signbit(fx_a) == signbit(fx_b) ){
        printf("Function has same sign in limits: f(%g) = %g f(%g) = %g.\n",x_a,fx_a,x_b,fx_b);
        assert(signbit(fx_a) != signbit(fx_b));
    }

    if( fabs(x_a - x_b) < epsilon ){
        printf("Limits are too close: |%g - %g| = %g < %g.\n",x_a,x_b,fabs(x_a - x_b),epsilon);
        assert(fabs(x_a - x_b) >= epsilon);
    }

    x_mean_previous = x_a;
    for(int i=1 ; i<=max_interations ; i++){
        x_mean = x_a + (x_b - x_a)/2;
        convergence = calculate_convergence(x_mean, x_mean_previous, relative_convergence);
        x_mean_previous = x_mean;
        
        if( convergence < convergence_tol){
            if( debug){
                fx_mean = func(x_mean);
                printf("Reached convergence after %i interations: f(%g) = %g.\n", i, x_mean, fx_mean);
            }
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
        printf("Convergence was not reached after %i interations: f(%g) = %g\n", max_interations, x_mean, fx_mean);
    return x_mean;
}

double find_root_bissection(double (*func)(double), double x_a, double x_b){
    return find_root_bissection_debug(func, x_a, x_b, default_convergence, true, default_interations, false);
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
    return max(0, round(n+0.5));
}

void test_bissection(double (*func)(double), double x_root, double x_a, double x_b, double convergence_tol, bool relative_convergence, int max_interations, bool debug, char *message){
    printf("\nTest Bissection Method: %s\n", message);
    int interations = estimate_bissection_interations(x_a, x_b, x_root, convergence_tol, false);
    printf("Estimated number of interations: %i\n", interations);
    double x_root_bissection = find_root_bissection_debug(func, x_a, x_b, convergence_tol, relative_convergence, max_interations, debug);
    print_error(" => Root", x_root, x_root_bissection);
}

// Root at x=0.3
double f_linear(double x){
    return -x*3 + 0.9;
}

// Root at x=0.3
double f_linear2(double x){
    return -x*3E5 + 0.9E5;
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

// Root at x= +-2
double f_exponential(double x){
    return exp(x*x-4) - 1;
}

int tests_bissection(){
    bool relative_convergence = false;
    int max_interations = 50;
    bool debug = true;

    test_bissection(f_linear, 0.3, 0, 2, 0.001, relative_convergence, max_interations, debug, "Linear function");
    test_bissection(f_linear, 0.3, 2, 0, 0.001, relative_convergence, max_interations, debug, "Linear function with inverted limits");
    test_bissection(f_linear, 0.3, 0, 2, 0.001, true, max_interations, debug, "Linear function with relative convergence");
    test_bissection(f_linear, 0.3, 0, 1, 0.001, relative_convergence, 5, debug, "Linear function with insufficient interations");
    test_bissection(f_linear, 0.3, 0., 0.4, 0.001, relative_convergence, max_interations, debug, "Linear function with early exit");
    test_bissection(f_linear, 0.3, 0.3, 1, 0.001, relative_convergence, max_interations, debug, "Linear function with root in x_a");
    test_bissection(f_linear, 0.3, 0, 0.3, 0.001, relative_convergence, max_interations, debug, "Linear function with root in x_b");
    test_bissection(f_linear, 0.3, 0.3-epsilon*100, 0.3+epsilon*100, 0.001, relative_convergence, max_interations, debug, "Linear function with domain close to root");
    test_bissection(f_linear, 0.3, 0.3-epsilon/3, 0.3+epsilon/3, 0.001, relative_convergence, max_interations, debug, "Linear function with domain very close to root");
    // test_bissection(f_linear, 0.3, 1, 2.3, 0.001, relative_convergence, max_interations, debug, "Linear function with error in [x_a,x_b]");
    // test_bissection(f_linear, 0.3, -2, 0., 0.001, relative_convergence, max_interations, debug, "Linear function with error in [x_a,x_b]");
    
    // test_bissection(f_linear2, 0.3, 0.3-epsilon/3, 0.3+epsilon/3, 0.001, relative_convergence, max_interations, debug, "Linear function #2 with 'small' domain");
    
    test_bissection(f_quadratic, -0.1, -0.25, 1, 0.001, relative_convergence, max_interations, debug, "Quadratic function, root#1");
    test_bissection(f_quadratic, -0.5, -0.25, -1, 0.001, relative_convergence, max_interations, debug, "Quadratic function, root#2");
    test_bissection(f_cos, pi/2, 0, 2, 0.001, relative_convergence, max_interations, debug, "Cossine function");
    test_bissection(f_trigonometric, 3./4*pi, 0, 5, 0.001, relative_convergence, max_interations, debug, "Trigonometric function");
    test_bissection(f_exponential, 2., 0, 10, 0.001, relative_convergence, max_interations, debug, "Exponential function");
}


// Minimum Curvature Method
double alfa(double theta1, double phi1, double theta2, double phi2){
    double x, y;
    x = sin( (theta2 - theta1)/2 );
    x *= x;
    y = sin( (phi2 - phi1)/2 );
    y *= y;
    y *= sin(theta1)*sin(theta2);
    x += y;
    x = sqrt(x);
    return 2* asin(x);
}

double f_alfa(double a){
    if( a<0.02){
        double x;
        double a2 = a*a;
        x = 1 + 32*a2/18;
        x = 1 + a2/168*x;
        x = 1+ a2/10*x;
        return 1+ a2/12*x;
    } else{
        return 2/a * tan(a/2);
    }
}

double deltaN_with_deltaSfa(double deltaSfa, double theta1, double phi1, double theta2, double phi2){
    return deltaSfa/2*(sin(theta1)*cos(phi1) + sin(theta2)*cos(phi2));
}

double deltaN_with_fa(double deltaS, double fa, double theta1, double phi1, double theta2, double phi2){
    return deltaN_with_deltaSfa(deltaS*fa, theta1, phi1, theta2, phi2);
}

double deltaN(double deltaS, double theta1, double phi1, double theta2, double phi2){
    double fa = f_alfa( alfa(theta1, phi1, theta2, phi2) );
    return deltaN_with_fa(deltaS, fa, theta1, phi1, theta2, phi2);
}

double deltaE_with_deltaSfa(double deltaSfa, double theta1, double phi1, double theta2, double phi2){
    return deltaSfa/2*(sin(theta1)*sin(phi1) + sin(theta2)*sin(phi2));
}

double deltaE_with_fa(double deltaS, double fa, double theta1, double phi1, double theta2, double phi2){
    return deltaE_with_deltaSfa(deltaS*fa, theta1, phi1, theta2, phi2);
}

double deltaE(double deltaS, double theta1, double phi1, double theta2, double phi2){
    double fa = f_alfa( alfa(theta1, phi1, theta2, phi2) );
    return deltaE_with_fa(deltaS, fa, theta1, phi1, theta2, phi2);
}


double deltaV_with_deltaSfa(double deltaSfa, double theta1, double theta2){
    return deltaSfa/2*(cos(theta1) + cos(theta2));
}

double deltaV_with_fa(double deltaS, double fa, double theta1, double theta2){
    return deltaV_with_deltaSfa(deltaS*fa, theta1, theta2);
}

double deltaV(double deltaS, double theta1, double phi1, double theta2, double phi2){
    double fa = f_alfa( alfa(theta1, phi1, theta2, phi2) );
    return deltaV_with_fa(deltaS, fa, theta1, theta2);
}

void tests_minimum_curvature(){
    double theta1, phi1, theta2, phi2;
    double a, fa;
    double dS, dN, dE, dV;

    printf("\n#### Minimum Curvatura Method Tests ####\n");

    printf("\nVertical well\n");
    theta1=0;
    phi1=0.88*pi;
    theta2=0;
    phi2=0.124*pi;
    dS=10;
    a = alfa(theta1, phi1, theta2, phi2);
    printf("  Alfa: calculated=%g  true=%g\n",a,0.);
    fa = f_alfa(a);
    printf("  f(alfa): calculated=%g  true=%g\n",fa,1.);
    dN = deltaN(dS, theta1, phi1, theta2, phi2);
    printf("  Delta N: calculated=%g  true=%g\n",dN,0.);
    dE = deltaE(dS, theta1, phi1, theta2, phi2);
    printf("  Delta E: calculated=%g  true=%g\n",dE,0.);
    dV = deltaV(dS, theta1, phi1, theta2, phi2);
    printf("  Delta V: calculated=%g  true=%g\n",dV,dS);

    printf("\nSlant well\n");
    theta1=pi/4;
    phi1= 0;

}

int main(){
    tests_bissection();
    tests_minimum_curvature();
    return 0;
}