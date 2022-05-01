#ifndef __IOUTIL_H__
#define __IOUTIL_H__

#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<iostream>
#include "Point.h"

using namespace std;

class IOUtil{
	
public:
    IOUtil(){}

    static void read_input_points(const char* fname, size_t& dim, vector<Point>& dataP);

    static void write_output_points(const char* fname, size_t dim, vector<Point>& dataP);

    static void readUtils(const char* fname, size_t dim, vector<Point>& Utils, size_t size);
   
    static void readAllUtils(const char* fname, size_t dim, vector<Point>& Utils);

    static void readRandomNumofUtils(const char* fname, size_t dim, vector<Point>& Utils, size_t size, int randomNumber);

    static void readRandomNumofUtilsNoD(const char* fname, size_t dim, vector<Point>& Utils, size_t size, int randomNumber);

    static void readFixedSizePoints(const char* fname, size_t& dim, vector<Point>& dataP, vector<int>& toInsert, size_t &size, vector<bool> &isDeleted);

    
};
#endif
