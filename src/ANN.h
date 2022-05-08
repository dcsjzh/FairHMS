#ifndef __ANN_H__
#define __ANN_H__

#include "Point.h"
#include <ANN/ANN.h>
#include <vector>

using namespace std;

class ANN{
  ANNbd_tree* tree;
  ANNpointArray pa;
  size_t data_size;
  size_t dim;
  size_t bucket_size;
  ANNsplitRule split;
  ANNshrinkRule shrink;

  private:
    void clean();
    void init();

  public:
    ANN();
    ~ANN();

    void insertPts(const vector<Point>& dataP);
    void getANNs(const vector<Point>& queryP, double epsilon, vector<size_t>& idxs);
    void getANNks(const vector<Point>& queryP, double epsilon, vector<size_t>& idxs, int r);
    void getkANNs(const vector<Point>& queryP, double epsilon, vector<vector<size_t>>& idxs, 
                                    vector<vector<double>> &k_dist,int r);
    void getSIngleANN(Point& queryP, double epsilon, size_t& idxs, double &dist);
    void getSInglekANN(Point& queryP, double epsilon, size_t& idxs, double &dist, int k);
    void getRANNs(const vector<Point> dataP, const vector<Point> Utils, Point& queryP, double epsilon, vector<size_t>& idxs, int k);
};

#endif
