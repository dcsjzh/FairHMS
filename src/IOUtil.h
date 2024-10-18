#ifndef __IOUTIL_H__
#define __IOUTIL_H__

#include <cstdlib>
#include <vector>
#include <iostream>
#include <set>
#include <iomanip>
#include <fstream>
#include "Point.h"
#include "RMSUtils.h"
using namespace std;

class IOUtil
{
public:
    IOUtil() {}
    static void read_input_points(const char *fname, size_t &dim, vector<Point> &dataP);
    static void read_input_functions(const char *fname, size_t &dim, vector<UtilityFunction> &FC);
    static void write_output_points(const char *fname, size_t dim, vector<Point> &dataP);
};
#endif
