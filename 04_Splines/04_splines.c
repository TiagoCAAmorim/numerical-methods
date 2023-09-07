/*
    Implementation of the bissection method to find root of 1D functions    
*/

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

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
    m[n] = 0.;
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
    double coef[4];
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
    double coef[4];
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
    double coef[4];
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

bool print_spline_table_file(struct tspline spline, double npoints, char *fname){
    double x0, x1, x;
    double y, yp, ypp, Sy;
    int j=1;

    FILE* file = fopen(fname, "w");

    if (file == NULL) {
        perror("Error creating file!");
        return false;
    }

    fprintf(file, "%16s\t%16s\t%16s\t%16s\t%16s\n","x","y","dy/dx","d2y/dx2","Int(y)");

    x0 = spline.x[0];
    x1 = spline.x[spline.npoints-1];
    for( int i = 0; i < npoints; i++) {
        x = x0 + (x1-x0) * i / (npoints - 1);
        if( x >=spline.x[j] && j < spline.npoints-1){
            x = spline.x[j];
            j++;
            if( x >spline.x[j]){
                i--;
            }
        };
        y = evaluate_spline(spline, x);
        yp = evaluate_spline_prime(spline, x);
        ypp = evaluate_spline_prime2(spline, x);
        Sy = evaluate_spline_integral(spline, x0, x);
        fprintf(file, "%16.10g\t%16.10g\t%16.10g\t%16.10g\t%16.10g\n", x, y, yp, ypp, Sy);
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
    assert(print_spline_table_file(spline, points_file, "Ex1.txt"));  
    
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
    assert(print_spline_table_file(spline, points_file, "Ex2.txt"));

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
    assert(print_spline_table_file(spline, points_file, "Ex3.txt"));
    
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
    assert(print_spline_table_file(spline, points_file, "Ex4.txt"));
}

int main(){
    test_splines();

    return 0;
}