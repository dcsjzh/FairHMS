#ifndef HS_H
#define HS_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "Point.h"
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include "IOUtil.h"
#include "RandUtil.h"
#include "MVE.h"
#include "ANN.h"
#include <boost/numeric/ublas/matrix.hpp>
#include "RMSUtils.h"
#include <ctime>
#include <queue> 
using namespace std;

void runHS(vector<Point> dataP, int r, int k,vector<Point> &result, vector<Point> curSky, double &time);

#endif