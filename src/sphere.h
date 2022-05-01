#ifndef SPHERE_H
#define SPHERE_H

#include "data_utility.h"
#include "search.h"
#include "operation.h"
#include"RMSUtils.h"

#include "lp.h"
#include"Point.h"

// The complete Sphere algorithm
point_set_t* sphereWSImpLP(point_set_t* point_set, int k);
void runSphere(vector<Point> dataP, int r, int k, vector<Point> &result, vector<Point> curSky, double &time);

#endif