#include "RMSUtils.h"
#include "RandUtil.h"
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <set>
#include "MVE.h"
#include <list>
#include <queue> 
#include <utility>

//initialize some gloabals
const vector<Point> * RMSUtils::ref_to_pts = NULL;
const Point * RMSUtils::ref_to_dir = NULL;
//Point temp_P;

size_t RMSUtils::dimension = 1;
boost::numeric::ublas::matrix<double> RMSUtils::transformation_matrix = boost::numeric::ublas::zero_matrix<double>(1,1);
boost::numeric::ublas::vector<double> RMSUtils::center = boost::numeric::ublas::zero_vector<double>(1);

/* returns when its good to use straightforward iteration
 * compared to make_heap
 */
bool smallk(size_t k, size_t n)
{
  assert(k <= n);

  //for very small n, make_heap is as fast
  if(n <= 1024)
    return false;

  return k <= ceil(log2((double) n));   
}

size_t RMSUtils::log_net_size(double sphere_radius, double net_radius, size_t dim)
{
  assert(sphere_radius > 0 && net_radius > 0);
  //cout << "get_net_size: sphere_radius = " << sphere_radius << ", net_radius = " << net_radius << ", dim = " << dim << endl;

  return 1 + (dim - 1) * (2  +  ceil(log2(sphere_radius / net_radius))) ;
}


size_t RMSUtils::log_random_net_size(double sphere_radius, double net_radius, double delta, size_t dim)
{
  
  assert(dim > 0 && sphere_radius > 0 && net_radius > 0 && 0 < delta &&
         delta < 1);

  size_t logM = RMSUtils::log_net_size(sphere_radius, net_radius/2, dim);

  //calculate log to base 2 of 2^{log M} * log_e( 2^{log M} / delta)
  //which is at most log M + log_2 ( logM - delta)
  double val = (double) logM  + log2 ((double) logM - log (delta));

  return ceil(val);
}

void RMSUtils::get_random_sphere_points(double sphere_radius, size_t dim, size_t N, 
                              vector<Point>& randomP,
                              bool first_orthant)
{

  Point p(dim);
  for(size_t i = 0; i < N; i++)
  {
    RandUtil::get_random_direction(dim,p);
    if(first_orthant)
      p = Point::abs(p);
    p.scale_to_length(sphere_radius);
    randomP.push_back(p);
  }
}



void RMSUtils::Max_Avg_Regret(double sphere_radius, size_t dim, size_t N, vector<Point>& dataP,
			      vector<Point>& R, double & MaxR, double & AvgR, double & perc80, int k){
  MaxR=0.0;
  double sum1=0.0;
  double maxQ1=-1.0;
  double maxPk=0.0;
  double temp;
  double maxRegret=0.0;
  double regret;
  double delta;
  //Point p(dim);
  
  vector<Point> random_dirs;
  RMSUtils::get_random_sphere_points(1.0, dim, N, random_dirs,true);
  vector<double> Scores;
  for(size_t i = 0; i < N; i++){   
    maxQ1=-1.0;
    for(size_t j=0;j<R.size();j++){
      temp=R[j].dotP(random_dirs[i]);
      if(temp>maxQ1){
	maxQ1=temp;
      }
    }
    size_t id;
    vector<size_t> usedPoints;
   
    for(size_t l=0;l<k;l++){
      maxPk=0.0;
      for(size_t j=0;j<dataP.size();j++){
	if(find(usedPoints.begin(), usedPoints.end(), j) != usedPoints.end()){
	  continue;
	}else{
	  temp=dataP[j].dotP(random_dirs[i]);	  
	  if(temp>maxPk){
	    maxPk=temp;
	    id=j;
	  }
	}
      }
      usedPoints.push_back(id);
    }
    if(maxPk>maxQ1){
      delta=maxQ1/maxPk;
      regret=1-delta;
      Scores.push_back(regret);
      sum1=sum1+regret;
      if(regret>maxRegret){
	maxRegret=regret;
      }	
    }else{
      Scores.push_back((double)0);
    }
  }
  if(Scores.size()>0){
    sort(Scores.begin(), Scores.end());
    int g=(int)(0.8*Scores.size()+0.5);
    perc80=Scores[g];
  }else{
    perc80=0;
  }
  MaxR=maxRegret;
  AvgR=sum1/N;
}




void RMSUtils::fast_Max_Avg_Regret(double sphere_radius, size_t dim, size_t N, vector<Point>& dataP,
			      vector<Point>& R, double & MaxR, double & AvgR, double & perc80, int k){
  MaxR=0.0;
  double sum1=0.0;
  double maxQ1=-1.0;
  double maxPk=0.0;
  double temp;
  double maxRegret=0.0;
  double regret;
  double delta;
  //Point p(dim);
  
  
  vector<Point> random_dirs;
  RMSUtils::get_random_sphere_points(1.0, dim, N, random_dirs,true);
  vector<double> Scores;
  for(size_t i = 0; i < N; i++){   
  	std::priority_queue<double> q;
    maxQ1=-1.0;
    for(size_t j=0;j<R.size();j++){
      temp=R[j].dotP(random_dirs[i]);
      if(temp>maxQ1){
	maxQ1=temp;
      }
    }
    size_t id;
    
    
    for(size_t j=0;j<dataP.size();j++){
    	q.push(dataP[j].dotP(random_dirs[i]));
	}
	
    for(int j=1;j<=k-1;j++){
    	q.pop();
	}    
    maxPk=q.top();
    
    if(maxPk>maxQ1){
      delta=maxQ1/maxPk;
      regret=1-delta;
      Scores.push_back(regret);
      sum1=sum1+regret;
      if(regret>maxRegret){
	maxRegret=regret;
      }	
    }else{
      Scores.push_back((double)0);
    }
    q=priority_queue<double>();
  }
  if(Scores.size()>0){
    sort(Scores.begin(), Scores.end());
    int g=(int)(0.8*Scores.size()+0.5);
    perc80=Scores[g];
  }else{
    perc80=0;
  }
  MaxR=maxRegret;
  AvgR=sum1/N;
}


void RMSUtils::Max_Avg_Regret2(double sphere_radius, size_t dim, size_t N, vector<Point>& dataP,
			      vector<Point>& R, double & MaxR, double & AvgR, vector<double>& perc, int k){
  MaxR=0.0;
  double sum1=0.0;
  double maxQ1=-1.0;
  double maxPk=0.0;
  double temp;
  double maxRegret=0.0;
  double regret;
  double delta;
  //Point p(dim);
  int cc=0;
  int ALL=0;
  
  vector<Point> random_dirs;
  RMSUtils::get_random_sphere_points(1.0, dim, N, random_dirs,true);
  vector<double> Scores;
  for(size_t i = 0; i < N; i++){   
    maxQ1=-1.0;
    for(size_t j=0;j<R.size();j++){
      temp=R[j].dotP(random_dirs[i]);
      if(temp>maxQ1){
	maxQ1=temp;
      }
    }
    size_t id;
    vector<size_t> usedPoints;
   
    for(size_t l=0;l<k;l++){
      maxPk=0.0;
      for(size_t j=0;j<dataP.size();j++){
	if(find(usedPoints.begin(), usedPoints.end(), j) != usedPoints.end()){
	  continue;
	}else{
	  temp=dataP[j].dotP(random_dirs[i]);	  
	  if(temp>maxPk){
	    maxPk=temp;
	    id=j;
	  }
	}
      }
      usedPoints.push_back(id);
    }
    if(maxPk>maxQ1){
      delta=maxQ1/maxPk;
      regret=1-delta;
      Scores.push_back(regret);
      sum1=sum1+regret;
      cc++;
      ALL++;
      if(regret>maxRegret){
	maxRegret=regret;
      }	
    }else{
      Scores.push_back((double)0);
      ALL++;
    }
  }
  sort(Scores.begin(), Scores.end());
  int g;
  for(int i=0;i<9;i++){
      g=(int)(0.1*(i+1)*Scores.size()+0.5);
      perc[i]=Scores[g];
  }
  perc[9]=Scores[Scores.size()-1];
  g=(int)(0.85*Scores.size()+0.5);
    perc[10]=Scores[g];
    g=(int)(0.95*Scores.size()+0.5);
    perc[11]=Scores[g];
  cout<<"cc is: "<<cc<<" out of "<<ALL<<"\n";
  MaxR=maxRegret;
  AvgR=sum1/N;
}

void RMSUtils::fast_Max_Avg_Regret2(double sphere_radius, size_t dim, size_t N, vector<Point>& dataP,
			      vector<Point>& R, double & MaxR, double & AvgR, vector<double>& perc, int k){
  MaxR=0.0;
  double sum1=0.0;
  double maxQ1=-1.0;
  double maxPk=0.0;
  double temp;
  double maxRegret=0.0;
  double regret;
  double delta;
  //Point p(dim);
  int cc=0;
  int ALL=0;
  
  vector<Point> random_dirs;
  RMSUtils::get_random_sphere_points(1.0, dim, N, random_dirs,true);
  vector<double> Scores;
  for(size_t i = 0; i < N; i++){ 
      std::priority_queue<double> q;
    maxQ1=-1.0;
    for(size_t j=0;j<R.size();j++){
      temp=R[j].dotP(random_dirs[i]);
      if(temp>maxQ1){
	maxQ1=temp;
      }
    }
    size_t id;

    for(size_t j=0;j<dataP.size();j++){
    	q.push(dataP[j].dotP(random_dirs[i]));
	}
	
    for(int j=1;j<=k-1;j++){
    	q.pop();
	}    
    maxPk=q.top();
    
    if(maxPk>maxQ1){
      delta=maxQ1/maxPk;
      regret=1-delta;
      Scores.push_back(regret);
      sum1=sum1+regret;
      cc++;
      ALL++;
      if(regret>maxRegret){
	maxRegret=regret;
      }	
    }else{
      Scores.push_back((double)0);
      ALL++;
    }
    q=priority_queue<double>();
  }
  sort(Scores.begin(), Scores.end());
  int g;
  for(int i=0;i<9;i++){
      g=(int)(0.1*(i+1)*Scores.size()+0.5);
      perc[i]=Scores[g];
  }
  perc[9]=Scores[Scores.size()-1];
  cout<<"cc is: "<<cc<<" out of "<<ALL<<"\n";
  MaxR=maxRegret;
  AvgR=sum1/N;
}

size_t RMSUtils::ndir_for_validation(size_t dim)
{
  size_t ndir;

  //number of random directions
  if(dim == 1)
    ndir = 2;
  else if(dim <= 3)
    ndir = 20000;
  else if(dim <= 4)
    ndir = 60000;
  else if(dim <= 6)
    ndir = 80000;//ndir = 32000;
  else if(dim <= 8)
    ndir = 120000;
  else if(dim <= 12){
    ndir = 170536;
    ndir = 200000;
  }else if(dim <= 16)
    ndir = 200000;
  else if(dim <= 20)
    ndir = 262144;
  else //dimension too large!
    ndir = 270500;

  return ndir;
}


/* returns true if idx1 is to be considered less than idx2 */
bool RMSUtils::heap_comp(size_t idx1, size_t idx2)
{
  assert(idx1 >= 0 && idx1 < ref_to_pts->size());
  assert(idx2 >= 0 && idx2 < ref_to_pts->size());

  if(idx1 == idx2)
    return false;

  /* evaluate the dot products with dir */
  double val1 = (*ref_to_pts)[idx1].dotP(*ref_to_dir);
  double val2 = (*ref_to_pts)[idx2].dotP(*ref_to_dir);

  if(val1 < val2)
    return true;

  if(val1 > val2)
    return false;

  //val1 == val2
  //if idx1 is later in the array then consider it smaller
  return idx1 > idx2;
}

void naive_rank_selection_dotp(const vector<Point>& pts,
                               const Point& dir,
                               size_t k,
                               vector<size_t>& topkI,
                               vector<double>& topkV)
{
  set < size_t > used;
  for(size_t i = 0; i < k; i++)
  {
    bool init = false;
    size_t midx;
    double maxv;
    for(size_t j = 0; j < pts.size(); j++)
    {
      //an already used index
      if(used.find(j) != used.end())
        continue;

      if(!init)
      {
        midx = j;
        maxv = pts[j].dotP(dir);
        init = true;
      }
      else
      {
        double curr = pts[j].dotP(dir);
        if(curr > maxv)
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
  
  for(size_t i = 0; i < k-1; i++)
    assert(topkV[i] >= topkV[i+1]);
  
  return;
}

void RMSUtils::rank_selection_dotp(const vector<Point>& pts,
                                const Point& dir,
                                size_t k,
				vector<size_t>& topkI,
                                vector<double>& topkV)
{
  assert(pts.size() > 0);

  assert(1 <= k && k <= pts.size());

  assert(topkI.size() == 0 && topkV.size() == 0);

  if(smallk(k,pts.size()))
  {
    naive_rank_selection_dotp(pts,dir,k,topkI,topkV);
    return;
  }

  /* insert indexes into a max heap by value of dot product with dir
   * then extract max k times
   */
  
  vector<size_t> MHeap;
  for(size_t i = 0; i < pts.size(); i++)
    MHeap.push_back(i);

  //first make sure comparator function works correctly
  ref_to_pts = &pts;
  ref_to_dir = &dir;
  make_heap(MHeap.begin(), MHeap.end(), heap_comp);

  //now extract the top k values and send back
  for(size_t i = 1; i <= k; i++)
  {
    pop_heap(MHeap.begin(), MHeap.end(),heap_comp);
    size_t Midx = MHeap.back();
    MHeap.pop_back();
    topkI.push_back(Midx);
    topkV.push_back(pts[Midx].dotP(dir));
  }

  for(size_t i = 0; i < k-1; i++)
    assert(topkV[i] >= topkV[i+1]);
  return;
}

void RMSUtils::get_fat_pointset(const vector<Point>& dataP, vector<Point>& fatP,
double& inner_rad, double& outer_rad)
{  
   MVEUtil::GetNormalizedMVE(dataP, 1.0, fatP, outer_rad, inner_rad);
}


void RMSUtils::Stavros_transform(const vector<Point>& dataP, size_t dim, vector<Point>& fatP)
{
  vector < double > maxV;

  for(size_t i = 0; i < dim; i++)
  {
    //find the max value of coordinate in dimension i
    double maxv = 0;
    for(size_t j = 0; j < dataP.size(); j++)
    {
      double ptjv = dataP[j].get_coordinate(i);
      if(ptjv < 0)
      {
        cout << "j = " << j << " ptjv = " << ptjv << endl;
        assert(false);
      }
      if(ptjv > maxv)
        maxv = ptjv;
    }

    assert(maxv > 0);

    maxV.push_back(maxv);
  }

  //pack up all of this as points
  for(size_t j = 0; j < dataP.size(); j++)
  {
    vector < double > coords;
    for(size_t i = 0; i < dim; i++)
    {
      coords.push_back( dataP[j].get_coordinate(i) / maxV[i]);
    }

    fatP.push_back(Point(dim,coords));
  }


}

void RMSUtils::get_fat_pointset2(const vector<Point>& dataP, vector<Point>& fatP, double& inner_rad, double& outer_rad)
{
  if(dataP.size() == 0)
  {
    inner_rad = 0;
    outer_rad = 0;
    return;
  }

  size_t dim = dataP[0].get_dimension();

  Stavros_transform(dataP, dim, fatP);

  double c = (double) dim + sqrt((double) dim);

  /* construct a point with these coordinates */
  vector < double > V(dim,1.0/c);

  Point center(dim,V);

  /* subtract this from all points of Stavros transform */
  for(size_t j = 0; j < fatP.size(); j++)
    fatP[j] = fatP[j] - center;							

  /* now we can choose outer radius to be 1.0 */
  outer_rad = 1.0;
  inner_rad = 1.0 / c;

  return;
}

//transfer pointset from vecrot to struct array for regret ratio calculation
point_set_t *RMSUtils::pointSetTransf(vector<Point> fatP)
{
	int number_of_points = fatP.size();
	int dim = fatP[0].get_dimension();
	point_set_t *point_set = alloc_point_set(number_of_points);

	int count = 0;
	for (vector<Point>::iterator it = fatP.begin(); it != fatP.end(); it++)
	{
		point_t *p = alloc_point(dim, (*it).getId());
		for (int j = 0; j < dim; j++)
		{
			p->coord[j] = (*it).get_coordinate(j);
		}
		point_set->points[count] = p;
		count++;
	}
	return point_set;
}

void RMSUtils::pointSetTransf(vector<Point> &curR, point_set_t *result4c)
{
	int dim = result4c->points[0]->dim;
	for (size_t i = 0; i < result4c->numberOfPoints; i++)
	{
		vector<double> coord;
		for (size_t j = 0; j < dim; j++)
		{
			coord.push_back(result4c->points[i]->coord[j]);
		}
		Point p(dim, result4c->points[i]->id, coord);
		curR.push_back(p);
	}
}