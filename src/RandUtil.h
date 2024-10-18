#ifndef __RANDUTIL_H__
#define __RANDUTIL_H__

#include "Point.h"
using namespace std;

class RandUtil{
public:
    static void get_random_direction(size_t dim, Point& rp);  
    static double randUnif(double lo, double hi);
};

#endif
