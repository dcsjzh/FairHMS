#ifndef __RMSUTILS_H__
#define __RMSUTILS_H__

#include "Point.h"
#include"data_utility.h"
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

class RMSUtils {
  private:
    static const vector<Point> * ref_to_pts;
    static const Point * ref_to_dir;
    static bool heap_comp(size_t idx1, size_t idx2);

  public:
    
    static size_t dimension;
    static boost::numeric::ublas::matrix<double> transformation_matrix;
    static boost::numeric::ublas::vector<double> center;

    /** 
      get the logarithm (to base 2) 
      of the size of a r-net on the sphere of radius R: 
      this is the smallest
      number of points such that every point is at a distance of at most
      r from some point of the net.
      One can calculate the size of the net to be 2(4R/r)^{d-1} where
      d is the dimension
    **/
    static size_t log_net_size(double sphere_radius, double net_radius,
                               size_t dim);


    /** 
      How many random points to generate on
      the sphere of radius R
      so that with high probability at least 1- delta it will be a r-net
      for the sphere?

      To answer this question: We take a r/2 net on the sphere
      Now imagine the M balls of radius r/2 covering the sphere.
      We generate Mlog_e (M / delta) points on the sphere uniformly at random.
      One can then show that we hit every of the balls and this happens
      with probability at least 1 - delta (so with high probability).
      Then, the random points thus taken form a lambda net on the
      sphere.

      This function returns the log of that number.
    **/
    static size_t log_random_net_size(double sphere_radius, double net_radius, double delta, size_t dim);

    /** 
      Get N random points on the sphere of radius sphere_radius in
      dim dimensions

      If the first_orthant flag is set then points in first orthant only are generated
    **/
    static void get_random_sphere_points(double sphere_radius, size_t dim, 
                                         size_t N, vector<Point>& randomP,
                                         bool first_orthant);


static void Max_Avg_Regret(double sphere_radius, size_t dim, size_t N, vector<Point>& dataP,
			      vector<Point>& R, double & MaxR, double & AvgR, double & perc80, int k);

static void fast_Max_Avg_Regret(double sphere_radius, size_t dim, size_t N, vector<Point>& dataP,
			      vector<Point>& R, double & MaxR, double & AvgR, double & perc80, int k);

static void Max_Avg_Regret2(double sphere_radius, size_t dim, size_t N, vector<Point>& dataP,
			      vector<Point>& R, double & MaxR, double & AvgR, vector<double>& perc, int k);

static void fast_Max_Avg_Regret2(double sphere_radius, size_t dim, size_t N, vector<Point>& dataP,
			      vector<Point>& R, double & MaxR, double & AvgR, vector<double>& perc, int k);
    /** 
      To check if some linear function exceeds another over the sphere
      check these many random directions
    **/
    static size_t ndir_for_validation(size_t dim);


    /**
      pts specifies a set of pts in R^dim for some dimension
      pts must not be empty

      dir is a direction

      k is an integer must be between 1 and size of idxs

      topkI denotes the top k indexes within pts 
      topkV denotes the top k values of ranking by dotp of the subset
    **/      
    static void rank_selection_dotp(const vector<Point>& pts,
                                    const Point& dir,
                                    size_t k,
				    vector<size_t>& topkI,
                                    vector<double>& topkV);

    /** 
      This function returns a pointset, fatP, with each element of
      fatP a certain linear transform of the point set dataP. The
      main property of this is that the width of the fatP in every
      direction is about the same, upto a constant factor. This is
      established via the returned values. An inner_rad and outer_rad
      are returned and it is guaranteed that the convex hull of fatP
      lies inside the ball B(0, outer_rad) while B(0, inner_rad) is
      inside it. The outer_rad and inner_rad are different only by
      a certain constant factor, depending on the dimension.
    **/
    static void get_fat_pointset(const vector<Point>& dataP, 
		          vector<Point>& fatP,
                          double& inner_rad, 
			  double& outer_rad);


    /**
      Implement Stavros's simple transform for points in first orthant 
    **/

    static void Stavros_transform(const vector<Point>& dataP, size_t dim, vector<Point>& fatP);
    /**
      Implements the simpler algorithm by Stavros - assumes that
      origin is present in point set and all points are in
      first orthant as well as the point set is full dimensional
   **/
    static void get_fat_pointset2(const vector<Point>& dataP,
                                  vector<Point>& fatP,
                                  double& inner_rad,
                                  double& outer_rad);

    static point_set_t *pointSetTransf(vector<Point> fatP);

    static void pointSetTransf(vector<Point> &curR, point_set_t* result4c);
};


#endif
