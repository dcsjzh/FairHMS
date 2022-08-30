#ifndef MYUTILS_H
#define MYUTILS_H

#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include "lp.h"
#include "Point.h"
#include "RMSUtils.h"
using namespace std;

struct fairofGroup 
{
    int lc, uc, ki;
};

void getSkyline(vector<Point> &dataP, vector<Point> &dataSky);
void writeToFile(vector<Point> &dataP, int dim, vector<unordered_map<string, int>> &groupToInt, unordered_map<int, fairofGroup> &fairnessConstraint,
                 string datasetPath, vector<int> &result, int k, int groupID, string algName, double time, double &MHR);
void writeToFileELD(vector<Point> &dataP, int dim, vector<unordered_map<string, int>> &groupToInt, unordered_map<int, fairofGroup> &fairnessConstraint,
                    string datasetPath, vector<int> &result, int k, int groupID, string algName, double time, double &MHR, double epsilon, double lambda, double deltaC, int setDelta);
void generateGroups(unordered_map<int, vector<Point>> &groupedDataP,
                   int groupID, vector<Point> dataP);
void readUtilityFunctions(string fileName, vector<Point> &utiFunClass, int m);
void readDataPoints(string fileName, vector<Point> &dataP, int &dim, int &groupDim,
                   vector<unordered_map<string, int>> &categoryToCodes);
vector<string> split(string src, char c);
unordered_map<int, fairofGroup> generateFairCons(unordered_map<int, vector<Point>> groupedDataP, vector<Point> dataP, int k);

#endif