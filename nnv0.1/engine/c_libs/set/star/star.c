#include <stdio.h>
#include<stdlib.h>
#include <errno.h>
#include <cblas.h>
//#include <stdarg.h>

/*
    c-library for star set and its methods for reachability analysis

    author: Dung Tran
    date: 4/24/2019

 */



typedef struct Star{

    /*
        S = c + V*a
          c : the center vector
          V : basis matrix
          a : is the predicate variables satisfy C * a <= d

        In this implementation, c = V[0], a[0] = 1, i.e., we put center vector in to the basis matrix

        Reference: 1) Star-based Reachability Analysis of Deep Neural Network, Hoang-Dung Tran, Submitted to FM2019

                   2) Simulation-Equivalent Reachability of Large Linear Systems with Inputs, S.Bak, CAV2017

    */


    /* parameters of the star set */

    int dim;  // dimension of star set
    int nv;  // number of predicate variables
    double *V; // basis matrix of star set, V[0] = c
    double *lb; // lower bound of predicate variables
    double *ub; // uper bound of predicate variables
    double *C; // constraint matrix of the predicate
    double *d; // constraint vector of the predicate


    void (*print)(); // print star information
    int (*is_empty)(); // check emptyness of a star
    int (*is_intersect)(); // check intersect with HalfSpace


} Star;

extern Star star(double *lb, double *ub);

int main() {

    double lb[2] = {1, 1.5};
    double ub[1] = {2};

    Star S = star(lb, ub);

    int l = sizeof(lb)/sizeof(lb[0]);

    printf("Dimension of the star S is: %d\n", S.dim);
    printf("Lower bound vector S is: %lf\n", S.lb[1]);
    printf("Upper bound vector S is: %lf\n", S.ub[1]);
    printf("Length of lower bound vector is: %d\n", l);

    return 0;

}


extern Star star(double *lb, double *ub){
    /*
       star constructor
    */

    Star S;

    int n1 = sizeof(lb)/sizeof(lb[0]);
    int n2 = sizeof(ub)/sizeof(ub[0]);

    printf("n1 = %d\n", n1);
    printf("n2 = %d\n", n2);

    if (n1 != n2) {


        perror("Inconsistent dimension");

    }

    S.lb = lb;
    S.ub = ub;

    return S;

}
