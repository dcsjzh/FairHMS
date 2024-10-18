#ifndef __RMSUTILS_H__
#define __RMSUTILS_H__

#include "Point.h"
#include <vector>

struct UtilityFunction{
    Point direction;
    double fmax;
};

using namespace std;

class RMSUtils {
private:
    static const vector<Point> * ref_to_pts;
    static const Point * ref_to_dir;
    static bool heap_comp(size_t idx1, size_t idx2);

public:
    static size_t dimension;
    static size_t log_net_size(double sphere_radius, double net_radius, size_t dim);
    static size_t log_random_net_size(double sphere_radius, double net_radius, double delta, size_t dim);
    static void get_random_sphere_points(double sphere_radius, size_t dim, size_t N, vector<Point>& randomP,bool first_orthant);
    static void get_random_utility_functions(double sphere_radius, size_t dim, size_t N, vector<UtilityFunction>& FunctionClass,bool first_orthant);
    static size_t ndir_for_validation(size_t dim);
    static void rank_selection_dotp(const vector<Point>& pts,const Point& dir,size_t k,vector<size_t>& topkI,vector<double>& topkV);
};


#endif
