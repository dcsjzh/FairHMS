#ifndef GREEDY_H
#define GREEDY_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <list>
#include "RandUtil.h"
#include <vector>
#include <set>
#include <iostream>
#include "IOUtil.h"
#include "MVE.h"
#include <boost/numeric/ublas/matrix.hpp>
#include "RMSUtils.h"
#include <ctime>
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <string>
#include "gurobi_c++.h"
#include "MyUtils.h"
#include "BiGreedy.h"
using namespace std;


void runGreedy(vector<Point> dataP, int r, int k,vector<Point> &result, vector<Point> curSky, double &time);

void runMatroidGreedy(vector<Point> dataP, int r, int k, vector<Point> &R, vector<Point> curSky, double &time, 
                                      int groupID,
                                      unordered_map<int, vector<Point>> &groupedDataP, 
                                      unordered_map<int, fairofGroup> &fairnessConstraint);

#endif