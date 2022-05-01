#include "ANN.h"
#include <cassert>
#include <stdexcept>

ANN::ANN(){
  init();
}

ANN::~ANN(){
  clean();
}

void ANN::init()
{
  tree = NULL;
  pa = NULL;
  data_size = 0;
  dim = 0;
  bucket_size = 1;
  split = ANN_KD_SUGGEST;
  shrink = ANN_BD_SUGGEST;
}

void ANN::clean()
{
  if(tree)
    delete tree;
  if(pa)
    annDeallocPts(pa);
  init();
}

void ANN::insertPts(const vector<Point>& dataP)
{
  clean();

  data_size = dataP.size();

  if(data_size == 0)
    return;

  //otherwise get dimension from that of first point
  dim = dataP[0].get_dimension();

  pa = annAllocPts(data_size,dim);

  for(size_t i = 0; i < data_size; i++)
  {
    for(size_t d = 0; d < dim; d++){
      pa[i][d] = dataP[i].get_coordinate(d);
    }
  }
  
  //code to fix bucket size - if bucket size is too small we run
  //into too deep a tree which leads to a fault

  bucket_size = 1 + ceil( (double) data_size / 10.0 );
  //cout << "InsertPts: data_size = " << data_size << " bucket_size = " << bucket_size << endl;
  tree = new ANNbd_tree(pa,data_size,dim,bucket_size,split,shrink);
}

void ANN::getANNs(const vector<Point>& queryP, double epsilon, 
vector<size_t>& idxs)
{
  //if this is an unitialized ANN, this is an incorrect usage
  if(tree == NULL)
    throw std::runtime_error("No points in ANN data structure");

  if(queryP.size() == 0)
    return;

  assert(dim == queryP[0].get_dimension());

  //allocate an array to hold the query points and copy over
  ANNpointArray qa = annAllocPts(queryP.size(), dim);

  for(size_t i = 0; i < queryP.size(); i++)
  {
    for(size_t d = 0; d < dim; d++)
      qa[i][d] = queryP[i].get_coordinate(d);
  }

  //allocate variables for holding indexes and distances
  ANNidx idx;
  ANNdist dist;

  //do the queries
  for(size_t i = 0; i < queryP.size(); i++)
  {
    //return only 1-nearest neighbor
    tree->annkSearch(qa[i], 1, &idx, &dist, epsilon);

    //now store the index
    idxs.push_back((size_t) idx);
  }

  annDeallocPts(qa);

  return; 
}

void ANN::getSIngleANN(Point& queryP, double epsilon, size_t& idxs, double &newDist)
{
  ANNidx idx;
  ANNdist dist;
  ANNpoint p = annAllocPt(queryP.get_dimension());
  for(size_t d = 0; d < dim; d++)
      p[d] = queryP.get_coordinate(d);
  tree->annkSearch(p, 1, &idx, &dist, epsilon);
  idxs = (size_t) idx;
  newDist = (double) dist;
}

void ANN::getkANNs(const vector<Point>& queryP, double epsilon, 
vector<vector<size_t>>& idxs, vector<vector<double>> &k_dist, int r)
{
  //if this is an unitialized ANN, this is an incorrect usage
  if(tree == NULL)
    throw std::runtime_error("No points in ANN data structure");

  if(queryP.size() == 0)
    return;

  assert(dim == queryP[0].get_dimension());

  //allocate an array to hold the query points and copy over
  ANNpointArray qa = annAllocPts(queryP.size(), dim);

  for(size_t i = 0; i < queryP.size(); i++)
  {
    for(size_t d = 0; d < dim; d++)
      qa[i][d] = queryP[i].get_coordinate(d);
  }

  //allocate variables for holding indexes and distances
  ANNidxArray idx;
  idx = new ANNidx[r];
  ANNdistArray dist;
  dist = new ANNdist[r];

  //do the queries
  for(size_t i = 0; i < queryP.size(); i++)
  {
    //return k-nearest neighbor
    tree->annkSearch(qa[i], r, idx, dist, epsilon);

    //now store the index
    vector<size_t> temp;
    vector<double> tempDist;
    for(int j=0;j<r;j++){
        temp.push_back((size_t) idx[j]);
        tempDist.push_back(sqrt((double)dist[j]));
    }
    idxs.push_back(temp);
    k_dist.push_back(tempDist);
    
    
  }
  delete [] idx;
  delete [] dist;
  annDeallocPts(qa);

  return; 
}

void ANN::getANNks(const vector<Point>& queryP, double epsilon, 
vector<size_t>& idxs, int r)
{
  //if this is an unitialized ANN, this is an incorrect usage
  if(tree == NULL)
    throw std::runtime_error("No points in ANN data structure");

  if(queryP.size() == 0)
    return;

  assert(dim == queryP[0].get_dimension());

  //allocate an array to hold the query points and copy over
  ANNpointArray qa = annAllocPts(queryP.size(), dim);

  for(size_t i = 0; i < queryP.size(); i++)
  {
    for(size_t d = 0; d < dim; d++)
      qa[i][d] = queryP[i].get_coordinate(d);
  }

  //allocate variables for holding indexes and distances
  ANNidxArray idx;
  idx = new ANNidx[r];
  ANNdistArray dist;
  dist = new ANNdist[r];

  //do the queries
  for(size_t i = 0; i < queryP.size(); i++)
  {
    //return only 1-nearest neighbor
    tree->annkSearch(qa[i], r, idx, dist, epsilon);

    //now store the index
    for(int j=0;j<r;j++){
        idxs.push_back((size_t) idx[j]);
      //  cout<<idx[j]<<","<<std::flush;
    }
    //cout<<"\n"<<std::flush;
  }
  delete [] idx;
  delete [] dist;
  annDeallocPts(qa);

  return; 
}

void ANN::getSInglekANN(Point& queryP, double epsilon, size_t& idxs, double &dist, int k)
{
  //allocate variables for holding indexes and distances
  ANNidxArray idx;
  idx = new ANNidx[k];
  ANNdistArray dists;
  dists = new ANNdist[k];
  ANNpoint p = annAllocPt(queryP.get_dimension());
  for(size_t d = 0; d < dim; d++)
      p[d] = queryP.get_coordinate(d);
  tree->annkSearch(p, k, idx, dists, epsilon);
  idxs = (size_t) idx[k-1];
  dist = (double) dists[k-1];
}

void ANN::getRANNs(const vector<Point> dataP, const vector<Point> Utils, 
                            Point& queryP, double epsilon, vector<size_t>& idxs, int k)
{
  for (size_t i = 0; i < Utils.size(); i++)
  {
    ANNpoint p = annAllocPt(queryP.get_dimension());
    for(size_t d = 0; d < dim; d++)
        p[d] = Utils[i].get_coordinate(d);
    ANNidxArray idx;
    idx = new ANNidx[k];
    ANNdistArray dists;
    dists = new ANNdist[k];
    tree->annkSearch(p, k, idx, dists, epsilon);
    for (size_t j = 0; j < k; j++)
    {
      if(dataP[idx[j]].getId() == queryP.getId())
      {
        idxs.push_back(Utils[i].getId());
        break;
      }
    }
  }
}
