#ifndef MYUTILS_H
#define MYUTILS_H

#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include "Point.h"
#include "RMSUtils.h"
using namespace std;

struct fairinCate 
{
    int lc, uc, ki;
};

void getSkyline(vector<Point> &dataP, vector<Point> &dataSky);
void writeToFile(vector<Point> &dataP, int dim, vector<unordered_map<string, int>> &categoryToInt, unordered_map<int, fairinCate> &partitionMatroid,
                 string dataset, vector<int> &result, int k, int cateID, string algName, double time, double &MHR);
void writeToFile2D(vector<Point> &dataP, int dim, vector<unordered_map<string, int>> &categoryToCodes, unordered_map<int, fairinCate> &partitionMatroid,
                   string dataset, vector<int> &result, int k, int cateID, string algName, double time, double mhr);
void groupGen(unordered_map<int, vector<Point>> &categorizedDataP,
                   int categoryId, vector<Point> dataP);
void readUtilityFunctions(string fileName, vector<Point> &utility, int m);
void readDataPoint(string fileName, vector<Point> &dataP, int &dim, int &categoryDim,
                   vector<unordered_map<string, int>> &categoryToCodes);
vector<string> split(string src, char c);
unordered_map<int, fairinCate> constructFairCons(unordered_map<int, vector<Point>> categorizedDataP, vector<Point> dataP, int k);

#endif