#ifndef LP_H
#define LP_H

#include <glpk.h>
#include "data_utility.h"

#include "operation.h"

/*
 * The linear programs for computing MRR
 */
double worstDirection(point_set_t *s, point_t* pt, double* &v);
double worstDirection(int index, point_set_t *s, point_t* pt, double* &v);
double worstDirection(int index, point_set_t *s, point_t* pt, float* &v);
double determinant(int n, double** a);


/*
 * Compute the MRR of a given set of points
 */
double evaluateLP(point_set_t *p, point_set_t* S, int VERBOSE);

#endif