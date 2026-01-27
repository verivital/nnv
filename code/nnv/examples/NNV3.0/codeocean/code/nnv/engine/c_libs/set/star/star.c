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

void test_create_star();
void test_multiple_indirection();


int main() {


    //test_create_star();
    test_multiple_indirection();

    return 0;

}


void test_create_star() {

    double lb[3] = {1, 1.5, 2};
    double ub[1] = {2};

    int l = sizeof(lb)/sizeof(lb[0]);
    int l1 = sizeof(lb);
    int l2 = sizeof(lb[0]);
    int u1 = sizeof(ub);
    int u2 = sizeof(ub[0]);


    printf("l1 = %d\n", l1);
    printf("l2 = %d\n", l2);
    printf("u1 = %d\n", u1);
    printf("u2 = %d\n", u2);

    printf("Address of lb %p\n", &lb);


    Star S1;
    Star S2;

    S1.dim = l;
    S1.lb = lb;
    S1.ub = ub;

    double lb2[2] = {2, 1.5};
    double ub2[2] = {2, 2};

    S2.dim = sizeof(lb2)/sizeof(lb2[0]);
    S2.lb = lb2;
    S2.ub = ub2;

    printf("Dimension of S1 is: %d\n", S1.dim);
    printf("Lower bound vector S1 is: %lf\n", S1.lb[2]);
    printf("3 element of the lower bound vector S1 is: %lf\n", S1.lb[3]);
    printf("Upper bound vector S1 is: %lf\n", S1.ub[0]);

    int N=100;
    Star *RS[N];

    RS[1] = &S1;
    RS[2] = &S2;

    printf("The first star dimension is: %d\n",RS[1]->dim);
    printf("The first star lower bound vector is: %lf\n", RS[1]->lb[0]);

    lb[0] = 0;
    l = 0;
    printf("The first star dimension is: %d\n",RS[1]->dim);
    printf("The first star lower bound vector is: %lf\n", RS[1]->lb[0]);

    printf("************What we should do*********** pass by value or pass by reference\n");

}

void test_multiple_indirection() {

    int a = 3;
    int *b;
    b = &a;
    int **c = &b;
    int ***d = &c;

    printf("a = %d \n", a);
    printf("*b = %d \n", *b);
    printf("**c = %d \n", **c);
    printf("***d = %d \n", ***d);

}
