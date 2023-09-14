/*
    Implementation of Spline interpolation to estimate intemediate values in VFP tables  
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

const double epsilon = 1E-7;
const int points_file = 201;

struct tspline{
    int npoints;
    double x[100];
    double x_prime[100];
    double a[100];
    double b[100];
    double c[100];
    double d[100];
};

void copy_farray(int npoints, double *original, double *new){
    for(int i=0; i<npoints; i++){
        new[i] = original[i];
    }
}

struct tspline build_natural_spline(int npoints, double *x, double *y){
    struct tspline spline;
    double h[100];
    double r[100];
    double m[100];
    double z[100];
    double alpha[100];
    int i, n;

    assert(npoints > 2);
    spline.npoints = npoints;
    copy_farray(npoints, x, spline.x);
    copy_farray(npoints, y, spline.a);

    n = npoints - 1;
    for(i=0; i<=n-1; i++){
        assert(x[i+1] > x[i]);
        h[i] = x[i+1] - x[i];
    }

    for(i=1; i<=n-1; i++){
        alpha[i] = 3. / h[i] * (y[i+1] - y[i]) - 3. / h[i-1] * (y[i] - y[i-1]);
    }

    r[0] = 1.;
    m[0] = 0.;
    z[0] = 0.;

    for(i=1; i<=n-1; i++){
        r[i] = 2. * (x[i+1] - x[i-1]) - h[i-1] * m[i-1];
        assert(fabs(r[i]) > epsilon);
        m[i] = h[i] / r[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / r[i];  
    }

    r[n] = 1.;
    z[n] = 0.;
    spline.c[n] = 0.;

    for(i=n-1; i>=0; i--){
        spline.c[i] = z[i] - m[i] * spline.c[i+1];
        spline.b[i] = (y[i+1] - y[i]) / h[i] - h[i] * (spline.c[i+1] + 2*spline.c[i]) / 3;
        spline.d[i] = (spline.c[i+1] - spline.c[i]) / 3 / h[i];
    }

    assert(fabs(spline.c[0]) < epsilon);

    return spline;
}

struct tspline build_fixed_spline(int npoints, double *x, double *y, double *yp){
    struct tspline spline;
    double h[100];
    double r[100];
    double m[100];
    double z[100];
    double alpha[100];
    int i, n;

    assert(npoints > 2);
    spline.npoints = npoints;
    copy_farray(npoints, x, spline.x);
    copy_farray(npoints, y, spline.a);

    n = npoints - 1;
    for(i=0; i<=n-1; i++){
        assert(x[i+1] > x[i]);
        h[i] = x[i+1] - x[i];
    }

    alpha[0] = 3. / h[0] * (y[1] - y[0]) - 3. * yp[0];
    alpha[n] = 3. * yp[n] - 3. / h[n-1] * (y[n] - y[n-1]);
    for(i=1; i<=n-1; i++){
        alpha[i] = 3. / h[i] * (y[i+1] - y[i]) - 3. / h[i-1] * (y[i] - y[i-1]);
    }

    r[0] = 2*h[0];
    m[0] = 0.5;
    z[0] = alpha[0] / r[0];

    for(i=1; i<=n-1; i++){
        r[i] = 2. * (x[i+1] - x[i-1]) - h[i-1] * m[i-1];
        assert(fabs(r[i]) > epsilon);
        m[i] = h[i] / r[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / r[i];  
    }

    r[n] = h[n-1] * (2 - m[n-1]);
    z[n] = (alpha[n] - h[n-1] * z[n-1]) / r[n];
    spline.c[n] = z[n];

    for(i=n-1; i>=0; i--){
        spline.c[i] = z[i] - m[i] * spline.c[i+1];
        spline.b[i] = (y[i+1] - y[i]) / h[i] - h[i] * (spline.c[i+1] + 2*spline.c[i]) / 3;
        spline.d[i] = (spline.c[i+1] - spline.c[i]) / 3 / h[i];
    }

    assert(fabs(spline.b[0] - yp[0]) < epsilon);

    return spline;
}

int find_interval(int npoints, double *x_array, double x_evaluate){
    for( int i=1; i<npoints-1; i++){
        if( x_evaluate <= x_array[i]){
            return i-1;
        }
    }
    return npoints-2;
}

double evaluate_linear(struct tspline spline, double x){
    double y;
    double h1, h2, h;
    int j;

    j = find_interval(spline.npoints, spline.x, x);
    h1 = x - spline.x[j];
    h2 = spline.x[j+1] - x;
    h = spline.x[j+1] - spline.x[j];

    y = spline.a[j] * h2/h + spline.a[j+1] * h1/h;

    return y;
}

double evaluate_spline(struct tspline spline, double x){
    double y=0;
    double h;
    double dx=1;
    double coef[4];
    int j;

    j = find_interval(spline.npoints, spline.x, x);
    h = x - spline.x[j];

    y = spline.a[j];
    dx *= h;
    y += spline.b[j] * dx;
    dx *= h;
    y += spline.c[j] * dx;
    dx *= h;
    y += spline.d[j] * dx;
    
    return y;
}

double evaluate_spline_prime(struct tspline spline, double x){
    double y=0;
    double h;
    double dx=1;
    int j;

    j = find_interval(spline.npoints, spline.x, x);
    h = x - spline.x[j];

    y = spline.b[j];
    dx *= h;
    y += spline.c[j] * dx * 2.;
    dx *= h;
    y += spline.d[j] * dx * 3.;
    
    return y;    
}

double evaluate_spline_prime2(struct tspline spline, double x){
    double y=0;
    double h;
    double dx=1;
    int j;

    j = find_interval(spline.npoints, spline.x, x);
    h = x - spline.x[j];

    y += spline.c[j] * dx * 2.;
    dx *= h;
    y += spline.d[j] * dx * 6.;
    
    return y;        
}

double evaluate_spline_integral_point(struct tspline spline, int interval, double x){
    double y=0;
    double h;
    double dx=1;
    int j;

    j = interval;
    h = x - spline.x[j];

    dx *= h;
    y = spline.a[j] * dx;
    dx *= h;
    y += spline.b[j] * dx / 2.;
    dx *= h;
    y += spline.c[j] * dx / 3.;
    dx *= h;
    y += spline.d[j] * dx / 4.;
    
    return y;
}

double evaluate_spline_integral(struct tspline spline, double x0, double x1){
    int j0, j1;
    double y;

    j0 = find_interval(spline.npoints, spline.x, x0);
    j1 = find_interval(spline.npoints, spline.x, x1);

    y = -evaluate_spline_integral_point(spline, j0, x0);

    if( j1 != j0){
        for( int j=j0; j<j1; j++){
            y += evaluate_spline_integral_point(spline, j, spline.x[j+1]);
        }
    }

    y += evaluate_spline_integral_point(spline, j1, x1);

    return y;
}

bool compare_double(double x1, double x2){
    if( fabs(x1 - x2) > epsilon){
        printf("Values are different: %16.10g - %16.10g = %16.10g", x1, x2, x1-x2);
        return false;
    }
    return true;
}

bool compare_farray(int npoints, double *array1, double *array2){
    for(int i=0; i<npoints; i++){
        if( fabs(array1[i] - array2[i]) > epsilon){
            return false;
        };
    }
    return true;
}

void delta_farray(int npoints, double *array1, double *array2, double array_out[]){
    for(int i=0; i<npoints; i++){
        array_out[i] = array1[i] - array2[i];
    }
}

void print_farray(char *pre_message, int npoints, double *x, char *post_message){
    assert(npoints > 0);
    printf("%s{", pre_message);
    for(int i=0; i<npoints; i++){
        printf("%16.10g",x[i]);
        if( i<npoints-1){
            printf(", ");
        }
    }
    printf("}%s", post_message);
}

void print_farray_compare_error(int npoints, double *array1, double *array2, char *message){
    double delta[100];
    char *s;
    s = " ";
    if( npoints > 3){
        s = "\n ";
    }
    printf("%s:",message);
    print_farray(s, npoints, array1, " -");
    print_farray(s, npoints, array2, " =");
    delta_farray(npoints, array1, array2, delta);
    print_farray(s, npoints, delta, "\n");
}

bool compare_splines(struct tspline spline1, struct tspline spline2){
    if( spline1.npoints != spline2.npoints){
        printf("Number of points in splines is different: %i and %i",spline1.npoints, spline2.npoints);
        return false;
    }
    if( !compare_farray(spline1.npoints, spline1.x, spline2.x)){
        print_farray_compare_error(spline1.npoints, spline1.x, spline2.x, "X arrays are different");
        return false;
    }
    if( !compare_farray(spline1.npoints-1, spline1.a, spline2.a)){
        print_farray_compare_error(spline1.npoints-1, spline1.a, spline2.a, "A arrays are different");
        return false;
    }
    if( !compare_farray(spline1.npoints-1, spline1.b, spline2.b)){
        print_farray_compare_error(spline1.npoints-1, spline1.b, spline2.b, "B arrays are different");
        return false;
    }
    if( !compare_farray(spline1.npoints-1, spline1.c, spline2.c)){
        print_farray_compare_error(spline1.npoints-1, spline1.c, spline2.c, "C arrays are different");
        return false;
    }
    if( !compare_farray(spline1.npoints-1, spline1.d, spline2.d)){
        print_farray_compare_error(spline1.npoints-1, spline1.d, spline2.d, "D arrays are different");
        return false;
    }
    return true;
}

bool check_spline_values(struct tspline spline, double npoints, double *x, double *y){
    double y_calc;
    for( int i=0; i<npoints; i++){
        y_calc = evaluate_spline(spline, x[i]);
        if( fabs(y[i] - y_calc) > epsilon){
            printf("%ith test failed at %g: %g != %g.\n",i,x[i],y[i],y_calc);
            return false;
        }
    }
    return true;
}

bool print_xy_table_file(char *fname, double npoints, double *x, double *y){
    FILE* file = fopen(fname, "w");

    if (file == NULL) {
        perror("Error creating file!");
        return false;
    }

    fprintf(file, "%16s\t%16s\n","x","y");

    for( int i = 0; i < npoints; i++) {
        fprintf(file, "%16.10g\t%16.10g\n", x[i], y[i]);
    }

    fclose(file);
    return true;
}

bool print_spline_table_file(struct tspline spline, double npoints, char *fname, int nxTests, double *xTests){
    double x0, x1, x;
    double y, yp, ypp, Sy, ylin;
    int j=1;
    double *xTry;
    int *n;

    if( nxTests > 0){
        xTry = xTests;
        n = &nxTests;
    } else{
        xTry = spline.x;
        n = &spline.npoints;
    }

    FILE* file = fopen(fname, "w");

    if (file == NULL) {
        perror("Error creating file!");
        return false;
    }

    fprintf(file, "%16s\t%16s\t%16s\t%16s\t%16s\t%16s\n","x","y","dy/dx","d2y/dx2","Int(y)","ylinear");

    x0 = xTry[0]; 
    x1 = xTry[*n-1]; 
    for( int i = 0; i < npoints; i++) {
        x = x0 + (x1-x0) * i / (npoints - 1);
        if( x >= xTry[j] && j < *n-1){
            x = xTry[j];
            j++;
            if( x >xTry[j]){
                i--;
            }
        };
        y = evaluate_spline(spline, x);
        yp = evaluate_spline_prime(spline, x);
        ypp = evaluate_spline_prime2(spline, x);
        Sy = evaluate_spline_integral(spline, x0, x);
        ylin = evaluate_linear(spline, x);
        fprintf(file, "%16.10g\t%16.10g\t%16.10g\t%16.10g\t%16.10g\t%16.10g\n", x, y, yp, ypp, Sy, ylin);
    }

    fclose(file);
    return true;
}

void test_splines(){
    struct tspline spline;
    struct tspline spline_true;
    double x[100];
    double y[100];
    double y_prime[100];
    double x_test[100];
    double y_test[100];
    int npoints;
    char message[100];

    sprintf(message, "Example 1 - 3 points (natural)");
    npoints = 3;
    x[0] = 1;
    y[0] = 2;
    x[1] = 2;
    y[1] = 3;
    x[2] = 3;
    y[2] = 5;
    spline_true.npoints = npoints;
    copy_farray(npoints, x, spline_true.x);
    copy_farray(npoints, y, spline_true.a);
    spline_true.b[0] = 3./4.;
    spline_true.b[1] = 3./2.;
    spline_true.c[0] = 0.;
    spline_true.c[1] = 3./4.;
    spline_true.d[0] = 1./4.;
    spline_true.d[1] = -1./4.;
    spline = build_natural_spline(npoints, x, y);
    if( !compare_splines(spline, spline_true)){
        printf("%s: Error in building spline.\n", message);
        assert(false);
    } else{
        printf("%s: No error in building spline.\n", message);
    }
    if( !compare_double(evaluate_spline_prime2(spline, x[2]), 0.)){
        printf("%s: Error in calculating last 2nd derivative.\n", message);
        assert(false);
    } else{
        printf("%s: No error in calculating last 2nd derivative.\n", message);
    }  
    assert(print_spline_table_file(spline, points_file, "Ex1.txt", 0, spline.x));  
    
    sprintf(message, "Example 2 - Exponential function (natural)");
    npoints = 4;
    for( int i=0; i<npoints; i++){
        x[i] = (double)i;
        y[i] = exp(x[i]);
    }
    spline_true.npoints = npoints;
    copy_farray(npoints, x, spline_true.x);
    copy_farray(npoints, y, spline_true.a);
    spline_true.b[0] = exp(1)-1-1./15*(-exp(3)+6*exp(2)-9*exp(1)+4);
    spline_true.b[1] = exp(2)-exp(1)-1./15*(2*exp(3)+3*exp(2)-12*exp(1)+7);
    spline_true.b[2] = exp(3)-exp(2)-1./15*(8*exp(3)-18*exp(2)+12*exp(1)-2);
    spline_true.c[0] = 0.;
    spline_true.c[1] = 1./5*(-exp(3)+6*exp(2)-9*exp(1)+4);
    spline_true.c[2] = 1./5*(4*exp(3)-9*exp(2)+6*exp(1)-1);
    spline_true.d[0] = 1./15*(-exp(3)+6*exp(2)-9*exp(1)+4);
    spline_true.d[1] = 1./3*(exp(3)-3*exp(2)+3*exp(1)-1);
    spline_true.d[2] = 1./15*(-4*exp(3)+9*exp(2)-6*exp(1)+1);
    spline = build_natural_spline(npoints, x, y);
    if( !compare_splines(spline, spline_true)){
        printf("%s: Error in building spline.\n", message);
        assert(false);
    } else{
        printf("%s: No error in building spline.\n", message);
    }
    x_test[0] = -1.3;
    y_test[0] = spline_true.a[0] + spline_true.b[0] * (-1.3) + spline_true.c[0] *1.3*1.3 + spline_true.d[0] *(-1.3)*1.3*1.3;
    x_test[1] = 0.75;
    y_test[1] = spline_true.a[0] + spline_true.b[0] * 0.75 + spline_true.c[0] *0.75*0.75 + spline_true.d[0] *0.75*0.75*0.75;
    x_test[2] = 1.3;
    y_test[2] = spline_true.a[1] + spline_true.b[1] * 0.3 + spline_true.c[1] *0.3*0.3 + spline_true.d[1] *0.3*0.3*0.3;
    x_test[3] = 2.87;
    y_test[3] = spline_true.a[2] + spline_true.b[2] * 0.87 + spline_true.c[2] *0.87*0.87 + spline_true.d[2] *0.87*0.87*0.87;
    x_test[4] = 4.4;
    y_test[4] = spline_true.a[2] + spline_true.b[2] * 2.4 + spline_true.c[2] *2.4*2.4 + spline_true.d[2] *2.4*2.4*2.4;

    if( !check_spline_values(spline, 5, x_test, y_test)){
        printf("%s: Error in calculating values with spline.\n", message);
        assert(false);
    } else{
        printf("%s: No error in calculating values with spline.\n", message);
    }
    if( !compare_double(evaluate_spline_prime2(spline, x[3]), 0.)){
        printf("%s: Error in calculating last 2nd derivative.\n", message);
        assert(false);
    } else{
        printf("%s: No error in calculating last 2nd derivative.\n", message);
    } 
    if( !compare_double(evaluate_spline_integral(spline, x[0],x[3]), 19.55228649)){
        printf("%s: Error in calculating integral.\n", message);
        assert(false);
    } else{
        printf("%s: No error in calculating integral.\n", message);
    }
    assert(print_spline_table_file(spline, points_file, "Ex2.txt", 0, spline.x));

    sprintf(message, "Example 3 - 3 points (fixed)");
    npoints = 3;
    x[0] = 1;
    y[0] = 2;
    x[1] = 2;
    y[1] = 3;
    x[2] = 3;
    y[2] = 5;
    y_prime[0] = 2;
    y_prime[2] = 1;
    spline_true.npoints = npoints;
    copy_farray(npoints, x, spline_true.x);
    copy_farray(npoints, y, spline_true.a);
    spline_true.b[0] = 2.;
    spline_true.b[1] = 3./2.;
    spline_true.c[0] = -5./2.;
    spline_true.c[1] = 2.;
    spline_true.d[0] = 3./2.;
    spline_true.d[1] = -3./2.;
    spline = build_fixed_spline(npoints, x, y, y_prime);
    if( !compare_splines(spline, spline_true)){
        printf("%s: Error in building spline.\n", message);
        assert(false);
    } else{
        printf("%s: No error in building spline.\n", message);
    }
    if( !compare_double(evaluate_spline_prime(spline, x[2]), y_prime[2])){
        printf("%s: Error in calculating last 1st derivative.\n", message);
        assert(false);
    } else{
        printf("%s: No error in calculating last 1st derivative.\n", message);
    }    
    assert(print_spline_table_file(spline, points_file, "Ex3.txt", 0, spline.x));
    
    sprintf(message, "Example 4 - Exponential function (fixed)");
    npoints = 4;
    for( int i=0; i<npoints; i++){
        x[i] = (double)i;
        y[i] = exp(x[i]);
        y_prime[i] = y[i];
    }
    spline_true.npoints = npoints;
    copy_farray(npoints, x, spline_true.x);
    copy_farray(npoints, y, spline_true.a);
    spline_true.b[0] = 1.0;
    spline_true.b[1] = 2.710162988;
    spline_true.b[2] = 7.326516343;
    spline_true.c[0] = 1./15*( 2*exp(3) -12*exp(2) +42*exp(1) -59);
    spline_true.c[1] = 1./15*(-4*exp(3) +24*exp(2) -39*exp(1) +28);
    spline_true.c[2] = 1./15*(14*exp(3) -39*exp(2) +24*exp(1)  -8);
    spline_true.d[0] = 0.2735993315; 
    spline_true.d[1] = 0.6951307906; 
    spline_true.d[2] = 2.019091618; 
    spline = build_fixed_spline(npoints, x, y, y_prime);
    if( !compare_splines(spline, spline_true)){
        printf("%s: Error in building spline.\n", message);
        assert(false);
    } else{
        printf("%s: No error in building spline.\n", message);
    }
    x_test[0] = -1.3;
    y_test[0] = spline_true.a[0] + spline_true.b[0] * (-1.3) + spline_true.c[0] *1.3*1.3 + spline_true.d[0] *(-1.3)*1.3*1.3;
    x_test[1] = 0.75;
    y_test[1] = spline_true.a[0] + spline_true.b[0] * 0.75 + spline_true.c[0] *0.75*0.75 + spline_true.d[0] *0.75*0.75*0.75;
    x_test[2] = 1.3;
    y_test[2] = spline_true.a[1] + spline_true.b[1] * 0.3 + spline_true.c[1] *0.3*0.3 + spline_true.d[1] *0.3*0.3*0.3;
    x_test[3] = 2.87;
    y_test[3] = spline_true.a[2] + spline_true.b[2] * 0.87 + spline_true.c[2] *0.87*0.87 + spline_true.d[2] *0.87*0.87*0.87;
    x_test[4] = 4.4;
    y_test[4] = spline_true.a[2] + spline_true.b[2] * 2.4 + spline_true.c[2] *2.4*2.4 + spline_true.d[2] *2.4*2.4*2.4;

    if( !check_spline_values(spline, 5, x_test, y_test)){
        printf("%s: Error in calculating values with spline.\n", message);
        assert(false);
    } else{
        printf("%s: No error in calculating values with spline.\n", message);
    }
        if( !compare_double(evaluate_spline_prime(spline, x[3]), y_prime[3])){
        printf("%s: Error in calculating last 1st derivative.\n", message);
        assert(false);
    } else{
        printf("%s: No error in calculating last 1st derivative.\n", message);
    }
    if( !compare_double(evaluate_spline_integral(spline, x[0],x[3]), 19.05964498)){
        printf("%s: Error in calculating integral.\n", message);
        assert(false);
    } else{
        printf("%s: No error in calculating integral.\n", message);
    } 
    assert(print_spline_table_file(spline, points_file, "Ex4.txt", 0, spline.x));
}


// ###### VFP Tables ######
struct tVFP{
    int number;
    double depth;
    int nLIQ;
    int nGLR;
    int nWCUT;
    int nLFG;
    int nWHP;
    double LIQ[10];
    double GLR[10];
    double WCUT[10];
    double LFG[10];
    double WHP[10];
    double BHP[10][10][10][10][10];
};

void read_VFP_file(char *fname, struct tVFP *vfp){
    FILE* file = fopen(fname, "r");
    
    int nLIQ, nGLR, nWCUT, nLFG, nWHP;
    int iLIQ, iGLR, iWCUT, iLFG, iWHP, iBHP;
    enum tcurrent {none, PTUBE, DEPTH, LIQ, GLR, WCUT, LFG, WHP, BHP};
    enum tcurrent current;

    if (file == NULL) {
        perror("Error reading file!");
        assert(false);
    }

    printf("Reading VFP file: %s...\n", fname);

    iLIQ = 0;
    iGLR = 0;
    iWCUT = 0;
    iLFG = 0;
    iWHP = 0;
    
    nLIQ = 0;
    nGLR = 0;
    nWCUT = 0;
    nLFG = 0;
    nWHP = 0;
    iBHP = 0;

    current = none;
    char message[30];
    while (fscanf(file, "%[^\n ] ", message) != EOF) {
        if(      strcmp(message,"*PTUBE1") == 0){ current = PTUBE;}
        else if( strcmp(message,"*DEPTH") == 0 ){ current = DEPTH;}
        else if( strcmp(message,"*LIQ"  ) == 0 ){ current = LIQ;}
        else if( strcmp(message,"*GLR"  ) == 0 ){ current = GLR;}
        else if( strcmp(message,"*WCUT" ) == 0 ){ current = WCUT;}
        else if( strcmp(message,"*LFG"  ) == 0 ){ current = LFG;}
        else if( strcmp(message,"*WHP"  ) == 0 ){ current = WHP;}
        else if( strcmp(message,"*BHP"  ) == 0 ){ current = BHP;}
        else{
            switch (current)
            {
            case PTUBE:
                vfp->number = atoi(message);
                break;
            case DEPTH:
                vfp->depth = atof(message);
                break;
            case LIQ:
                vfp->LIQ[nLIQ] = atof(message);
                nLIQ += 1;
                break;           
            case GLR:
                vfp->GLR[nGLR] = atof(message);
                nGLR += 1;
                break;           
            case WCUT:
                vfp->WCUT[nWCUT] = atof(message);
                nWCUT += 1;
                break;           
            case LFG:
                vfp->LFG[nLFG] = atof(message);
                nLFG += 1;
                break;           
            case WHP:
                vfp->WHP[nWHP] = atof(message);
                nWHP += 1;
                break;           
            case BHP:
                if( iLIQ == 0){iLIQ = atoi(message);}
                else if( iGLR == 0){iGLR = atoi(message);}
                else if( iWCUT == 0){iWCUT = atoi(message);}
                else if( iLFG == 0){iLFG = atoi(message);}
                else{
                    vfp->BHP[iLIQ-1][iGLR-1][iWCUT-1][iLFG-1][iWHP] = atof(message);
                    iWHP += 1;
                    iBHP += 1;
                    if( iWHP == nWHP){
                        iLIQ=0;
                        iGLR=0;
                        iWCUT=0;
                        iLFG=0;
                        iWHP=0;
                    }
                }
                break;           
            default:
                assert(false);
                break;
            }
        }
    }
    printf("Table #%i at %g:\n", vfp->number, vfp->depth);
    printf(" %i LIQ:  [%10.4g;%10.4g]\n",nLIQ,  vfp->LIQ [0], vfp->LIQ [nLIQ-1]);
    printf(" %i GLR:  [%10.4g;%10.4g]\n",nGLR,  vfp->GLR [0], vfp->GLR [nGLR-1]);
    printf(" %i WCUT: [%10.4g;%10.4g]\n",nWCUT, vfp->WCUT[0], vfp->WCUT[nWCUT-1]);
    printf(" %i LFG:  [%10.4g;%10.4g]\n",nLFG,  vfp->LFG [0], vfp->LFG [nLFG-1]);
    printf(" %i WHP:  [%10.4g;%10.4g]\n",nWHP,  vfp->WHP [0], vfp->WHP [nWHP-1]);
    printf(" %i BHP\n",iBHP);
    fflush(stdout);
    assert( iBHP == (nLIQ*nGLR*nWCUT*nLFG*nWHP));
    vfp->nLIQ  = nLIQ;
    vfp->nGLR  = nGLR;
    vfp->nWCUT = nWCUT;
    vfp->nLFG  = nLFG;
    vfp->nWHP  = nWHP;
    fclose(file);
}

struct tPoints{
    double x[100];
    double y[100];
    int n;
};

struct tTestPoints{
    struct tPoints data;
    struct tPoints test;
};

void SeparatePoints(struct tPoints *points, struct tTestPoints *testPoints, int first, int step){
    int is_test = first;

    testPoints->data.n = 0;
    testPoints->test.n = 0;
    for( int i=0; i<points->n; i++){
        if( is_test == 0){
            testPoints->test.x[testPoints->test.n] = points->x[i];
            testPoints->test.y[testPoints->test.n] = points->y[i];
            testPoints->test.n += 1;
            is_test = step;
        } else{
            testPoints->data.x[testPoints->data.n] = points->x[i];
            testPoints->data.y[testPoints->data.n] = points->y[i];
            testPoints->data.n += 1;
        }
        is_test -= 1;
    }
}

void get_data_from_vfp(struct tVFP *vfp, int iGLR, int iWCUT, int iLFG, int iWHP, struct tPoints *points){
    points->n = vfp->nLIQ;
    for( int i=0; i<vfp->nLIQ; i++){
        points->x[i] = vfp->LIQ[i];
        points->y[i] = vfp->BHP[i][iGLR][iWCUT][iLFG][iWHP];
    }
}

double calculate_mean(const double data[], const int filter[], int n, int filterValue) {
    double sum = 0.0;
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (filter[i]==filterValue) {
            sum += data[i];
            count++;
        }
    }

    if (count == 0) {
        printf("No valid elements found for mean calculation.\n");
        return 0.0;
    }

    return sum / count;
}

double calculate_StdDev(const double data[], const int filter[], int n, int filterValue) {
    double mean = calculate_mean(data, filter, n, filterValue);
    double sum = 0.0;
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (filter[i] == filterValue) {
            double diff = data[i] - mean;
            sum += diff * diff;
            count++;
        }
    }

    if (count == 0) {
        printf("No valid elements found for standard deviation calculation.\n");
        return 0.0; 
    }

    return sqrt(sum / count);
}

void test_vfp(char *fname, struct tVFP *vfp, bool extrapolate, bool fixed, double beta0, double betan){
    struct tPoints points;
    struct tTestPoints testPoints;
    struct tspline spline;
    double yTrue;
    double ySpline;
    double yLinear;
    double fp[100];
    int n=0;
    int iLIQ;
    FILE* file = fopen(fname, "w");
    int iLIQstart, iLIQend;
    double eSpline[10000];
    double eLinear[10000];
    int LIQ[10000];
    double meanESpl, stdDevESpl;
    double meanELin, stdDevELin;

    if( extrapolate){
        iLIQstart = 0;
        iLIQend = vfp->nLIQ;
    } else{
        iLIQstart = 1;
        iLIQend = vfp->nLIQ - 1;
    }


    if (file == NULL) {
        perror("Error creating file!");
        printf("Couldn't write file %s\n",fname);
        assert(false);
    }
    fprintf(file, "%4s\t%4s\t%4s\t%4s\t%4s\t%-16s\t%-16s\t%-16s\t%-16s\t%-16s\n","LIQ","GLR","WCUT","LFG","WHP","True","Spline","eSpline","Linear","eLinear");

    for( int iGLR=0; iGLR<vfp->nGLR; iGLR++){
        for( int iWCUT=0; iWCUT<vfp->nWCUT; iWCUT++){
            for( int iLFG=0; iLFG<vfp->nLFG; iLFG++){
                for( int iWHP=0; iWHP<vfp->nWHP; iWHP++){
                    get_data_from_vfp(vfp, iGLR, iWCUT, iLFG, iWHP, &points);
                    
                    if( fixed){
                        fp[0] = (points.y[1] - points.y[0]) / (points.x[1] - points.x[0]) * beta0;
                        fp[vfp->nLIQ-2] = (points.y[vfp->nLIQ-1] - points.y[vfp->nLIQ-2]) / (points.x[vfp->nLIQ-1] - points.x[vfp->nLIQ-2]) * betan;
                    }

                    for( int j=iLIQstart; j<iLIQend; j++){
                        SeparatePoints(&points, &testPoints, j, 99);
                        if( fixed){
                            spline = build_fixed_spline(testPoints.data.n, testPoints.data.x, testPoints.data.y, fp);
                        } else{
                            spline = build_natural_spline(testPoints.data.n, testPoints.data.x, testPoints.data.y);
                        }
                        for( int i=0; i<testPoints.test.n; i++){
                            yTrue = testPoints.test.y[i];
                            ySpline = evaluate_spline(spline, testPoints.test.x[i]);
                            yLinear = evaluate_linear(spline, testPoints.test.x[i]);
                            iLIQ = j;
                            fprintf(file, "%4i\t%4i\t%4i\t%4i\t%4i\t",iLIQ+1,iGLR+1,iWCUT+1,iLFG+1,iWHP+1);
                            fprintf(file, "%16.10g\t%16.10g\t%16.10g\t%16.10g\t%16.10g\n", yTrue, ySpline, yTrue - ySpline, yLinear, yTrue - yLinear);

                            eSpline[n] = yTrue - ySpline; 
                            eLinear[n] = yTrue - yLinear; 
                            LIQ[n] = iLIQ; 
                            n += 1;
                        }
                    }
                }
            }
        }
    }
    fclose(file);


    printf("Error Statistics by LIQ Value\n");
    printf("%4s\t%16s\t%16s\t%16s\t%16s\n","LIQ","meanESpl", "stdDevESpl", "meanELin", "stdDevELin");
    for( int j=iLIQstart; j<iLIQend; j++){
        meanESpl = calculate_mean(eSpline, LIQ, n, j);
        meanELin = calculate_mean(eLinear, LIQ, n, j);
        stdDevESpl = calculate_StdDev(eSpline, LIQ, n, j);
        stdDevELin = calculate_StdDev(eLinear, LIQ, n, j);
        
        printf("%4i\t%16.10g\t%16.10g\t%16.10g\t%16.10g\n",j,meanESpl, stdDevESpl, meanELin, stdDevELin);
    }

}

void test_vfp_print_table(char *fname, struct tVFP *vfp, int iGLR, int iWCUT, int iLFG, int iWHP, int iLIQ){
    struct tPoints points;
    struct tTestPoints testPoints;
    struct tspline spline;
    
    get_data_from_vfp(vfp, iGLR, iWCUT, iLFG, iWHP, &points);
    SeparatePoints(&points, &testPoints, iLIQ, 99);
    spline = build_natural_spline(testPoints.data.n, testPoints.data.x, testPoints.data.y);
    if( !print_spline_table_file(spline, 100, fname, points.n, points.x)){
        printf("Couldn't create file %s.", fname);
    }
}

void test_vfp_print_table_true(char *fname, struct tVFP *vfp, int iGLR, int iWCUT, int iLFG, int iWHP){
    struct tPoints points;
    
    get_data_from_vfp(vfp, iGLR, iWCUT, iLFG, iWHP, &points);
    if( !print_xy_table_file(fname, points.n, points.x, points.y)){
        printf("Couldn't create file %s.", fname);
    }
}

void tests_vfp_interpolation(){
    struct tVFP vfp;
    char folder[] = "C:/Users/tiago.LENOVO-I7/Unicamp/2023.02/IM253_MetodosNumericos/Trabalhos/04_Splines/vfp/";
    char fname[256];

    sprintf(fname,"%s%s",folder,"P1.inc");
    read_VFP_file(fname, &vfp);
    sprintf(fname,"%s%s",folder,"P1.txt");
    test_vfp(fname, &vfp, true, false,1,1);
    sprintf(fname,"%s%s",folder,"P1_noExtrap.txt");
    test_vfp(fname, &vfp, false, false,1,1);

    sprintf(fname,"%s%s", folder, "P1_Ok.txt");
    test_vfp_print_table(fname, &vfp, 1-1, 1-1, 1-1, 3-1, 5-1);
    sprintf(fname,"%s%s",folder,"P1_Ok_true.txt");
    test_vfp_print_table_true(fname, &vfp, 1-1, 1-1, 1-1, 3-1);
    sprintf(fname,"%s%s", folder, "P1_BadSpline.txt");
    test_vfp_print_table(fname, &vfp, 2-1, 6-1, 1-1, 1-1, 3-1);
    sprintf(fname,"%s%s",folder,"P1_BadSpline_true.txt");
    test_vfp_print_table_true(fname, &vfp, 2-1, 6-1, 1-1, 1-1);
    sprintf(fname,"%s%s", folder, "P1_BadLinear.txt");
    test_vfp_print_table(fname, &vfp, 1-1, 1-1, 2-1, 1-1, 3-1);
    sprintf(fname,"%s%s",folder,"P1_BadLinear_true.txt");
    test_vfp_print_table_true(fname, &vfp, 1-1, 1-1, 2-1, 1-1);
    sprintf(fname,"%s%s", folder, "P1_Bad.txt");
    test_vfp_print_table(fname, &vfp, 2-1, 6-1, 1-1, 2-1, 2-1);
    sprintf(fname,"%s%s",folder,"P1_Bad_true.txt");
    test_vfp_print_table_true(fname, &vfp, 2-1, 6-1, 1-1, 2-1);
    sprintf(fname,"%s%s", folder, "P1_ExtrapLower.txt");
    test_vfp_print_table(fname, &vfp, 1-1, 4-1, 2-1, 3-1, 1-1);
    sprintf(fname,"%s%s",folder,"P1_ExtrapLower_true.txt");
    test_vfp_print_table_true(fname, &vfp, 1-1, 4-1, 2-1, 3-1);
    sprintf(fname,"%s%s", folder, "P1_ExtrapUpper.txt");
    test_vfp_print_table(fname, &vfp, 5-1, 4-1, 3-1, 2-1, 6-1);
    sprintf(fname,"%s%s",folder,"P1_ExtrapUpper_true.txt");
    test_vfp_print_table_true(fname, &vfp, 5-1, 4-1, 3-1, 2-1);


    sprintf(fname,"%s%s",folder,"P2.inc");
    read_VFP_file(fname, &vfp);
    sprintf(fname,"%s%s",folder,"P2.txt");
    test_vfp(fname, &vfp, true, false,1,1);

    sprintf(fname,"%s%s",folder,"P3.inc");
    read_VFP_file(fname, &vfp);
    sprintf(fname,"%s%s",folder,"P3.txt");
    test_vfp(fname, &vfp, true, false,1,1);

    sprintf(fname,"%s%s",folder,"P4.inc");
    read_VFP_file(fname, &vfp);
    sprintf(fname,"%s%s",folder,"P4.txt");
    test_vfp(fname, &vfp, true, false,1,1);

    sprintf(fname,"%s%s",folder,"P5.inc");
    read_VFP_file(fname, &vfp);
    sprintf(fname,"%s%s",folder,"P5.txt");
    test_vfp(fname, &vfp, true, false,1,1);

    sprintf(fname,"%s%s",folder,"P6.inc");
    read_VFP_file(fname, &vfp);
    sprintf(fname,"%s%s",folder,"P6.txt");
    test_vfp(fname, &vfp, true, false,1,1);

    sprintf(fname,"%s%s",folder,"P7.inc");
    read_VFP_file(fname, &vfp);
    sprintf(fname,"%s%s",folder,"P7.txt");
    test_vfp(fname, &vfp, true, false,1,1);

    sprintf(fname,"%s%s",folder,"P8.inc");
    read_VFP_file(fname, &vfp);
    sprintf(fname,"%s%s",folder,"P8.txt");
    test_vfp(fname, &vfp, true, false,1,1);

    sprintf(fname,"%s%s",folder,"P9.inc");
    read_VFP_file(fname, &vfp);
    sprintf(fname,"%s%s",folder,"P9.txt");
    test_vfp(fname, &vfp, true, false,1,1);

    sprintf(fname,"%s%s",folder,"P10.inc");
    read_VFP_file(fname, &vfp);
    sprintf(fname,"%s%s",folder,"P10.txt");
    test_vfp(fname, &vfp, true, false,1,1);


    sprintf(fname,"%s%s",folder,"P1.inc");
    read_VFP_file(fname, &vfp);
    
    double beta0;
    for( int i=0; i<35; i++){
        sprintf(fname,"%s%s%i%s",folder,"P1_fixed_beta0_",i,".txt");
        beta0 = 0.6 + 4.*i/40.;
        printf("Beta0 = %g\n", beta0);
        test_vfp(fname, &vfp, false, true,beta0,1.);
    }

    double betan;
    for( int i=0; i<10; i++){
        sprintf(fname,"%s%s%i%s",folder,"P1_fixed_betan_",i,".txt");
        betan = 0.4 + 1.*i/10.;
        printf("%i. Betan = %g\n", i, betan);
        test_vfp(fname, &vfp, false, true,1.8,betan);
    }

    sprintf(fname,"%s%s",folder,"P1_fixed.txt");
    test_vfp(fname, &vfp, false, true,1.8,0.9);

}

int main(){
    test_splines();
    tests_vfp_interpolation();
    return 0;
}