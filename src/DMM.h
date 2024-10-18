#ifndef DMM_H
#define DMM_H

#include "data_utility.h"
#include "operation.h"
#include <algorithm>
#include <set>
#include"Point.h"
#include"RMSUtils.h"

point_set_t* DMM(point_set_t* point_set, int k);

point_set_t* DMM_Greedy(point_set_t* point, int k);

void runDMMGreedy(vector<Point> dataP, int r, int k,vector<Point> &result, vector<Point> curSky, double &time);

void runDMMRRMS(vector<Point> dataP, int r, int k,vector<Point> &result, vector<Point> curSky, double &time);

#endif
