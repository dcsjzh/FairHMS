#include <cmath>
#include <cassert>
#include <algorithm>
#include <set>
#include <list>
#include <queue>
#include <utility>
#include "RMSUtils.h"
#include "RandUtil.h"

const vector<Point> *RMSUtils::ref_to_pts = NULL;
const Point *RMSUtils::ref_to_dir = NULL;

template <typename T>
T RandT(T _min, T _max)
{
    T temp;
    if (_min > _max)
    {
        temp = _max;
        _max = _min;
        _min = temp;
    }
    return rand() / (double)RAND_MAX * (_max - _min) + _min;
}

void RMSUtils::get_random_sphere_points(double sphere_radius, size_t dim, size_t N, vector<Point> &randomP, bool first_orthant)
{
    Point p(dim);
    for (size_t i = 0; i < N; i++)
    {
        RandUtil::get_random_direction(dim, p);
        if (first_orthant)
        {
            p = Point::abs(p);
        }
        p.scale_to_length(sphere_radius);
        randomP.push_back(p);
    }
}

void RMSUtils::get_random_utility_functions(double sphere_radius, size_t dim, size_t N, vector<UtilityFunction> &FunctionClass, bool first_orthant)
{
    Point p(dim);
    for (size_t i = 0; i < N; i++)
    {
        UtilityFunction utilityfunction;
        RandUtil::get_random_direction(dim, p);
        if (first_orthant)
        {
            p = Point::abs(p);
        }
        p.scale_to_length(sphere_radius);
        utilityfunction.direction = p;
        utilityfunction.fmax = 0;
        FunctionClass.push_back(utilityfunction);
    }
}

size_t RMSUtils::ndir_for_validation(size_t dim)
{
    size_t ndir;
    if (dim == 1)
        ndir = 2;
    else if (dim <= 3)
        ndir = 20000;
    else if (dim <= 4)
        ndir = 60000;
    else if (dim <= 6)
        ndir = 80000;
    else if (dim <= 8)
        ndir = 120000;
    else if (dim <= 12)
        ndir = 170536;
    else if (dim <= 16)
        ndir = 200000;
    else if (dim <= 20)
        ndir = 262144;
    else
        ndir = 270500;

    return ndir;
}

bool smallk(size_t k, size_t n)
{
    assert(k <= n);
    if (n <= 1024)
        return false;
    return k <= ceil(log2((double)n));
}

size_t RMSUtils::log_net_size(double sphere_radius, double net_radius, size_t dim)
{
    assert(sphere_radius > 0 && net_radius > 0);
    return 1 + (dim - 1) * (2 + ceil(log2(sphere_radius / net_radius)));
}

size_t RMSUtils::log_random_net_size(double sphere_radius, double net_radius, double delta, size_t dim)
{
    assert(dim > 0 && sphere_radius > 0 && net_radius > 0 && 0 < delta && delta < 1);
    size_t logM = RMSUtils::log_net_size(sphere_radius, net_radius / 2, dim);
    double val = (double)logM + log2((double)logM - log(delta));
    return ceil(val);
}

bool RMSUtils::heap_comp(size_t idx1, size_t idx2)
{
    assert(idx1 >= 0 && idx1 < ref_to_pts->size());
    assert(idx2 >= 0 && idx2 < ref_to_pts->size());
    if (idx1 == idx2)
        return false;
    double val1 = (*ref_to_pts)[idx1].dotP(*ref_to_dir);
    double val2 = (*ref_to_pts)[idx2].dotP(*ref_to_dir);
    if (val1 < val2)
        return true;
    if (val1 > val2)
        return false;
    return idx1 > idx2;
}

void naive_rank_selection_dotp(const vector<Point> &pts,const Point &dir,size_t k,vector<size_t> &topkI,vector<double> &topkV)
{
    set<size_t> used;
    for (size_t i = 0; i < k; i++)
    {
        bool init = false;
        size_t midx;
        double maxv;
        for (size_t j = 0; j < pts.size(); j++)
        {
            if (used.find(j) != used.end())
                continue;
            if (!init)
            {
                midx = j;
                maxv = pts[j].dotP(dir);
                init = true;
            }
            else
            {
                double curr = pts[j].dotP(dir);
                if (curr > maxv)
                {
                    midx = j;
                    maxv = curr;
                }
            }
        }
        topkI.push_back(midx);
        topkV.push_back(maxv);
        used.insert(midx);
    }
    for (size_t i = 0; i < k - 1; i++)
        assert(topkV[i] >= topkV[i + 1]);
    return;
}

void RMSUtils::rank_selection_dotp(const vector<Point> &pts,const Point &dir,size_t k,vector<size_t> &topkI,vector<double> &topkV)
{
    assert(pts.size() > 0);
    assert(1 <= k && k <= pts.size());
    assert(topkI.size() == 0 && topkV.size() == 0);
    if (smallk(k, pts.size()))
    {
        naive_rank_selection_dotp(pts, dir, k, topkI, topkV);
        return;
    }
    vector<size_t> MHeap;
    for (size_t i = 0; i < pts.size(); i++)
        MHeap.push_back(i);
    ref_to_pts = &pts;
    ref_to_dir = &dir;
    make_heap(MHeap.begin(), MHeap.end(), heap_comp);
    for (size_t i = 1; i <= k; i++)
    {
        pop_heap(MHeap.begin(), MHeap.end(), heap_comp);
        size_t Midx = MHeap.back();
        MHeap.pop_back();
        topkI.push_back(Midx);
        topkV.push_back(pts[Midx].dotP(dir));
    }
    for (size_t i = 0; i < k - 1; i++)
        assert(topkV[i] >= topkV[i + 1]);
    return;
}