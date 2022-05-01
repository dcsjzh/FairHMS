#ifndef GREEDYK_H
#define GREEDYK_H

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
using namespace std;


void runGreedyK(vector<Point> dataP, int r, int k,vector<Point> &result, vector<Point> curSky, double &time);
void runGreedy(vector<Point> dataP, int r, int k,vector<Point> &result, vector<Point> curSky, double &time);

#endif