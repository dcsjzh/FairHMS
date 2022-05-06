#ifndef DATA_UTILITY_H
#define DATA_UTILITY_H

#include<math.h>
#include<float.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include <vector>
#include <iostream>
#include <iterator>

#define COORD_TYPE			double
#define DIST_TYPE			double
#define PI					3.1415926
#define INF					100000000
#define	MAX_FILENAME_LENG	256
#define EQN_EPS				1e-9


// data structure for storing points.
typedef struct point
{
	int			dim;
	COORD_TYPE*	coord;

	int			id;

}	point_t;

// data structure for storing point set.
typedef struct point_set
{
	int numberOfPoints;
	point_t **points;
}	point_set_t;

// date structure for computing catersian product
typedef std::vector<double> Vi;
typedef std::vector<Vi> Vvi;

// point & point_set manipulation
point_t* alloc_point(int dim);
point_t* alloc_point(int dim, int id);
void release_point( point_t* &point_v);
point_set_t* alloc_point_set(int numberOfPoints);
void release_point_set(point_set_t* &point_set_v, bool clear);
void print_point(point_t* point_v);
void print_point_set(point_set_t* point_set_v);

#endif