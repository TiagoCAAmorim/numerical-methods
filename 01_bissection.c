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

const bool stop_on_error = true;

bool check_error(bool ok_condition, char *error_message){
    if( !ok_condition){
        printf(error_message);
        if( stop_on_error){
            assert(false);
        }
    }
    return ok_condition;
}

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

void print_error(char *message, double true_value, double calculated_value, double error_limit){
    double relative_error;
    relative_error = calculate_convergence(true_value, calculated_value, true);
    printf("%s: True=%g Calculated=%g Error=%g", message, true_value, calculated_value, relative_error);
    if( relative_error > error_limit){
        printf("  <= ######### Atention ###########");
    }
    printf("\n");
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
        if( debug){
            double x_trial;
            printf("Sample of function results in the provided domain:\n");
            for(int i=0; i<=20; i++){
                x_trial = x_a + (x_b - x_a)*i/20;
                printf("  f(%g) = %g\n", x_trial, func(x_trial));
            }
        }
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
            printf("%3i.  f(%-10g) = %-10g\tConvergence=%-10g\n",i,x_mean,fx_mean, convergence);
        
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
    print_error(" => Root", x_root, x_root_bissection, convergence_tol);
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
double angle_subtraction(double a2, double a1){
    double dif = a2 - a1;
    if( dif < -pi){
        dif = dif + 2*pi;
    } else if( dif > pi){
        dif = dif - 2*pi;
    }
    return dif;
}

double alfa(double theta1, double phi1, double theta2, double phi2){
    double x, y;
    x = sin( angle_subtraction(theta2, theta1)/2 );
    x *= x;
    y = sin( angle_subtraction(phi2, phi1)/2 );
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

struct tangle{
    double cos, sin, rad;
};

double calculate_rad(struct tangle angle, bool debug){
    check_error( fabs(angle.cos) <= 1., "## Failed |cos(angle)| <= 1.\n");
    check_error( fabs(angle.sin) <= 1., "## Failed |sin(angle)| <= 1.\n");
    
    angle.cos = min(1,max(angle.cos,-1));
    angle.sin = min(1,max(angle.sin,-1));
    
    double a_cos = acos(angle.cos);
    double a_sin = asin(angle.sin);

    if( angle.sin < 0){
        a_cos = 2*pi - a_cos;
        if( angle.cos > 0){
            a_sin = 2*pi + a_sin;
        }
    }
    if( angle.cos < 0){
        a_sin = pi - a_sin;
    }

    double a = (a_sin + a_cos)/2;
    if( debug){
        double convergence = calculate_convergence(a, a_cos, true);
        printf("Convergence error in angle calculation (%4g.pi rad): %g\n", a/pi, convergence);
    }
    return a;
}

struct tangle calculate_theta2(double deltaV, double deltaSfa, double cos_theta1){
    struct tangle theta2;
    theta2.cos = 2*deltaV / deltaSfa - cos_theta1;
    if( check_error( fabs(theta2.cos) <= 1., "## Failed |cos(theta2)| <= 1.\n")){
        theta2.sin = sqrt(1 - theta2.cos * theta2.cos);
    } else{
        theta2.sin = 0;
    }
    return theta2;    
}

struct tangle calculate_phi2_deltaE_zero(double deltaN, double sin_theta1, double cos_phi1, double sin_phi1, double sin_theta2){
    struct tangle phi2;
    if( fabs(sin_theta2)<epsilon && fabs(sin_theta1)<epsilon){
        phi2.sin = -sin_phi1;
    } else{
        if( check_error( fabs(sin_theta2) > epsilon, "## Failed sin(theta2) !=0.\n")){
            phi2.sin = - sin_theta1 / sin_theta2 * sin_phi1;
        } else{
            phi2.sin = -sin_phi1;
        }
    }
    if( check_error( fabs(phi2.sin) <= 1, "## Failed |sin(phi2)| <= 1.\n")){
        phi2.cos = sqrt(1 - phi2.sin*phi2.sin);
        if( cos_phi1 < 0){
            phi2.cos = -phi2.cos;
        }
    } else{
        phi2.cos = 0;
    }
    return phi2;
}

struct tangle calculate_phi2_deltaN_zero(double deltaE, double sin_theta1, double cos_phi1, double sin_phi1, double sin_theta2){
    struct tangle phi2;
    if( fabs(sin_theta2)<epsilon && fabs(sin_theta1)<epsilon){
        phi2.cos = -cos_phi1;
    } else{
        check_error( fabs(sin_theta2) > epsilon, "## Failed sin(theta2) != 0.\n");
        phi2.cos = - sin_theta1 / sin_theta2 * cos_phi1;
    }
    phi2.sin = sqrt(1 - phi2.sin*phi2.sin);
    if( sin_phi1 > 0){
        phi2.sin = -phi2.sin;
    }
    return phi2;
}

struct tangle calculate_phi2(double deltaE, double deltaN, double sin_theta1, double cos_phi1, double sin_phi1, double sin_theta2){
    struct tangle phi2;

    if( fabs(sin_theta2) < epsilon){
        phi2.cos = cos_phi1;
        phi2.sin = sin_phi1;
    } else if( fabs(deltaE) < epsilon && fabs(deltaN) < epsilon){
        phi2.cos = -cos_phi1;
        phi2.sin = -sin_phi1;
    } else if( fabs(deltaE) < epsilon){
        phi2 = calculate_phi2_deltaE_zero(deltaN, sin_theta1, cos_phi1, sin_phi1, sin_theta2);
    } else if( fabs(deltaN) < epsilon){
        phi2 = calculate_phi2_deltaN_zero(deltaE, sin_theta1, cos_phi1, sin_phi1, sin_theta2);
     } else{
        double deltaEpsilon = deltaN * sin_phi1 - deltaE * cos_phi1;
        double deltaEpsilon_sin_theta1 = deltaEpsilon * sin_theta1;
        double deltaH2 = deltaE*deltaE + deltaN*deltaN;
        double deltaBeta2 = deltaH2 * sin_theta2*sin_theta2 - deltaEpsilon_sin_theta1*deltaEpsilon_sin_theta1;
        
        check_error( deltaBeta2 >= 0, "## Failed deltaBeta2 >= 0.\n");
        deltaBeta2 = max(0, deltaBeta2);
        double deltaBeta = sqrt(deltaBeta2);

        double deltaH2_sin_theta2 = deltaH2 * sin_theta2;
        if( !check_error( fabs(deltaH2_sin_theta2) > epsilon, "## Failed deltaH2_sin_theta2 != 0.\n")){
            deltaH2_sin_theta2 = deltaH2*0.001;
        }

        phi2.sin = (-deltaN * deltaEpsilon_sin_theta1 + deltaE*deltaBeta) / deltaH2_sin_theta2;
        phi2.cos = ( deltaE * deltaEpsilon_sin_theta1 + deltaN*deltaBeta) / deltaH2_sin_theta2;
    }
    return phi2;
}

double calculate_deltaSfa(double deltaE, double deltaN, double deltaV, double cos_theta1, double sin_theta1, double cos_phi1, double sin_phi1, double cos_theta2, double sin_theta2, double cos_phi2, double sin_phi2){
    double aE = sin_theta1 * sin_phi1 + sin_theta2 * sin_phi2;
    double aN = sin_theta1 * cos_phi1 + sin_theta2 * cos_phi2;
    double aV = cos_theta1 + cos_theta2;
    return 2 * sqrt( (deltaE*deltaE + deltaN*deltaN + deltaV*deltaV) / (aE*aE + aN*aN + aV*aV) );
}

double calculate_deltaSfa_aproximate(double deltaE, double deltaN, double deltaV, double cos_theta1, double sin_theta1, double cos_phi1, double sin_phi1, double deltaSfa){
    struct tangle theta2 = calculate_theta2(deltaV, deltaSfa, cos_theta1);
    struct tangle phi2 = calculate_phi2(deltaE, deltaN, sin_theta1, cos_phi1, sin_phi1, theta2.sin);
    return calculate_deltaSfa(deltaE, deltaN, deltaV, cos_theta1, sin_theta1, cos_phi1, sin_phi1, theta2.cos, theta2.sin, phi2.cos, phi2.sin);
}

double calculate_deltaSfa_error(double deltaE, double deltaN, double deltaV, double cos_theta1, double sin_theta1, double cos_phi1, double sin_phi1, double deltaSfa){
    double deltaSfa_calc = calculate_deltaSfa_aproximate( deltaE, deltaN, deltaV, cos_theta1, sin_theta1, cos_phi1, sin_phi1, deltaSfa);
    return deltaSfa_calc - deltaSfa;
}

double dE, dN, dV, cos_theta1, sin_theta1, cos_phi1, sin_phi1;
void define_well_path_data(double deltaE, double deltaN, double deltaV, double theta1, double phi1){
    dE = deltaE;
    dN = deltaN;
    dV = deltaV;
    dE = deltaE;
    cos_theta1 = cos(theta1);
    sin_theta1 = sin(theta1);
    cos_phi1 = cos(phi1);
    sin_phi1 = sin(phi1);
}
double calculate_defined_deltaSfa_error(double deltaSfa){
    return calculate_deltaSfa_error(dE, dN, dV, cos_theta1, sin_theta1, cos_phi1, sin_phi1, deltaSfa);
}

void test_MCM_formulas(char *message, double deltaS, double theta1, double phi1, double theta2, double phi2, double true_alfa, double true_deltaE, double true_deltaN, double true_deltaV, bool report_alfa_dEdNdV){

    printf("\n%s\n", message);

    double a = alfa(theta1, phi1, theta2, phi2);
    double dE = deltaE(deltaS, theta1, phi1, theta2, phi2);
    double dN = deltaN(deltaS, theta1, phi1, theta2, phi2);
    double dV = deltaV(deltaS, theta1, phi1, theta2, phi2);

    if( report_alfa_dEdNdV){
        print_error("  Alfa",true_alfa,a, 0.01);
        print_error("  Delta E", true_deltaE,dE, 0.01);
        print_error("  Delta N", true_deltaN,dN, 0.01);
        print_error("  Delta V", true_deltaV,dV, 0.01);
    } else{
        printf("  Calculated displacement: dE=%g dN=%g dV=%g\n", dE, dN, dV);
        printf("  Calculated alfa=%g\n",a);
        printf("  Calculated f(alfa)=%g\n",f_alfa(a));
    }

    double deltaSfa = deltaS * f_alfa(a);
    struct tangle theta2_ = calculate_theta2(dV, deltaSfa, cos(theta1));
    print_error("  theta2", theta2, calculate_rad(theta2_, false), 0.01);

    struct tangle phi2_ = calculate_phi2(dE, dN, sin(theta1), cos(phi1), sin(phi1), sin(theta2));
    print_error("  phi2", phi2, calculate_rad(phi2_, false), 0.01);

    double deltaSfa_calc = calculate_deltaSfa(dE, dN, dV, cos(theta1), sin(theta1), cos(phi1), sin(phi1), cos(theta2), sin(theta2), cos(phi2), sin(phi2));
    double dS = deltaSfa_calc / f_alfa(a);
    print_error("  Delta S", deltaS, dS, 0.01);

    printf("Calculate theta2, phi2 and dS from dE, dN, dV, theta1 and phi1.\n");
    define_well_path_data( dE, dN, dV, theta1, phi1);
    double deltaSfa_min = sqrt(dE*dE + dN*dN + dV*dV);
    double deltaSfa_max = deltaSfa_min * f_alfa(0.95*pi);
    if( dV > 0){
        deltaSfa_min = max(deltaSfa_min, 2*dV/(cos_theta1+1) );
    } else if( dV < 0){
        deltaSfa_min = max(deltaSfa_min, 2*dV/(cos_theta1-1) );
    }

    deltaSfa_calc = find_root_bissection_debug(calculate_defined_deltaSfa_error, deltaSfa_min, deltaSfa_max, 0.001, true, 100, true);
        
    print_error("  Delta S x f(alfa)", deltaSfa, deltaSfa_calc, 0.01);
    theta2_ = calculate_theta2(dV, deltaSfa_calc, cos(theta1));
    print_error("  theta2", theta2, calculate_rad(theta2_, false), 0.01);
    phi2_ = calculate_phi2(dE, dN, sin(theta1), cos(phi1), sin(phi1), theta2_.sin);
    print_error("  phi2", phi2, calculate_rad(phi2_, false), 0.01);
    a = alfa(theta1, phi1, calculate_rad(theta2_, false), calculate_rad(phi2_, false));
    dS = deltaSfa_calc / f_alfa(a);
    print_error("  Delta S", deltaS, dS, 0.01);
}

void tests_minimum_curvature(){
    double theta1, phi1, theta2, phi2;
    double a;
    double dS, dE, dN, dV;
    struct tangle theta2_, phi2_;
    double theta2_calc, phi2_calc;
    double dSfa, dS_calc;

    printf("\n#### Minimum Curvatura Method Tests ####\n");

    printf("\nAngle calculation\n");
    for(int i=0; i<=10; i++){
        a = 0.2*i*pi;
        theta2_.cos = cos(a);
        theta2_.sin = sin(a);
        theta2_.rad = calculate_rad(theta2_, true);
        print_error(" Angle calculation",a,theta2_.rad, 0.01);
    }

    dS=10;
    theta1=0;    phi1=0.88*pi;
    theta2=0;    phi2=0.88*pi;
    a=0.;
    dE=0., dN=0., dV=dS;
    test_MCM_formulas("Vertical well", dS, theta1, phi1, theta2, phi2, a, dE, dN, dV, true);

    dS=10;
    theta1=pi/4;    phi1=pi/6;
    theta2=pi/4;    phi2=pi/6;
    a=0.;
    dE=dS*sin(pi/4)*sin(pi/6);
    dN=dS*sin(pi/4)*cos(pi/6);
    dV=dS*cos(pi/4);
    test_MCM_formulas("Slant straight well", dS, theta1, phi1, theta2, phi2, a, dE, dN, dV, true);
    
    dS=10*pi/2;
    theta1=pi/2;    phi1=3*pi/2;
    theta2=pi/2;    phi2=0.;
    a=pi/2;
    dE=-10;    dN=10;    dV=0.;
    test_MCM_formulas("1/4 circle horizontal well", dS, theta1, phi1, theta2, phi2, a, dE, dN, dV, true);

    dS=10*pi/2;
    theta1=pi/2;    phi1=pi/4;
    theta2=pi/2;    phi2=7*pi/4.;
    a=pi/2;
    dE=0.;    dN=10*sqrt(2);    dV=0.;
    test_MCM_formulas("1/4 circle deltaE=0 horizontal well", dS, theta1, phi1, theta2, phi2, a, dE, dN, dV, true);
    
    dS=10*pi/2;
    theta1=0;        phi1=0.1871*pi;
    theta2=pi/2;     phi2=0.;
    a=pi/2;
    dE=0.;    dN=10.;    dV=10.;
    test_MCM_formulas("1/4 circle alligned north vertical to horizontal well", dS, theta1, phi1, theta2, phi2, a, dE, dN, dV, true);
    
    dS=10*pi/2;
    theta1=pi/6;         phi1=0.;
    theta2=theta1+pi/2;  phi2=0.;
    a=pi/2;
    dE=0.;    dN=10.*(cos(pi/6)+sin(pi/6));    dV=10.*(cos(pi/6)-sin(pi/6));
    test_MCM_formulas("1/4 circle alligned north well 'going up'", dS, theta1, phi1, theta2, phi2, a, dE, dN, dV, true);
    
    dS=10*pi/2;
    theta1=pi/4;         phi1=0.;
    theta2=theta1+pi/2;  phi2=0.;
    a=pi/2;
    dE=0.;    dN=10.*(cos(theta1)+sin(theta1));    dV=10.*(cos(theta1)-sin(theta1));
    test_MCM_formulas("1/4 circle alligned north well 'going up' #2", dS, theta1, phi1, theta2, phi2, a, dE, dN, dV, true);
    
    dS=10*pi/2;
    theta1=pi/3;         phi1=0.;
    theta2=theta1+pi/2;  phi2=0.;
    a=pi/2;
    dE=0.;    dN=10.*(cos(theta1)+sin(theta1));    dV=10.*(cos(theta1)-sin(theta1));
    test_MCM_formulas("1/4 circle alligned north well 'going up' #3", dS, theta1, phi1, theta2, phi2, a, dE, dN, dV, true);

    dS=10*pi/2;
    theta1=0;        phi1=0.1871*pi;
    theta2=pi/2;     phi2=pi/6;
    a=pi/2;
    dE=10.*sin(pi/6);    dN=10.*cos(pi/6);    dV=10.;
    test_MCM_formulas("1/4 circle 30o north vertical to horizontal well", dS, theta1, phi1, theta2, phi2, a, dE, dN, dV, true);
    
    dS=10.;
    theta1=pi/10;      phi1=pi/4;
    theta2=pi/6;       phi2=7*pi/4;
    a=0.;
    dE=0.;    dN=0.;    dV=0.;
    test_MCM_formulas("3D path well", dS, theta1, phi1, theta2, phi2, a, dE, dN, dV, false);   
}

int main(){
    // tests_bissection();
    tests_minimum_curvature();
    return 0;
}