#include"hs.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include"HittingSet.h"
#include "RandUtil.h"
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include "IOUtil.h"
#include "MVE.h"
#include <boost/numeric/ublas/matrix.hpp>
#include "RMSUtils.h"
#include <queue> 
#include <ctime>

using namespace std;

float clockMhs=(float)0.0;
float clockM2hs=(float)0.0;

bool validate_hitting_set(const vector < set < int > >& Sets, const vector < int>& HS)
{
 
  for(size_t i = 0; i < Sets.size(); i++)
  {
    bool hit = false;
    for(size_t j = 0; j < HS.size(); j++)
    {
      set < int > :: iterator itA;
      itA = Sets[i].find(HS[j]);
      if(itA != Sets[i].end())
      {
        hit = true;
	break;
      }
    }

    if(!hit)
    {
      return false;
    }
  }

  return true;
}

void test_examples()
{
  set < int > A, B, C, D;
  vector < set < int > > Sets1;
  vector < int > HS1;

  A.insert(1);
  A.insert(2);
  A.insert(3);

  B.insert(2);

  C.insert(2);
  C.insert(3);

  D.insert(1);
  D.insert(3);
  D.insert(4);

  Sets1.push_back(A);
  Sets1.push_back(B);
  Sets1.push_back(C);
  Sets1.push_back(D);

  HSApprox::get_hs_approximation<int>(Sets1, HS1);

  bool valid = validate_hitting_set(Sets1,HS1);
  if(valid)
  {
    cout << "Hitting set for example 1 is valid " << endl;
    for(size_t i = 0; i < HS1.size(); i++)
      cout << HS1[i] << "  ";
    cout << endl;
  }
  else
  {
    cout << "Hitting set for example 1 is not valid " << endl;
    exit(1);
  }

  //a more substantial example
  vector < set < int > > Sets2;
  for(size_t i = 0; i < 100000; i++)
  {
    set < int > S;

    for(int j = 1; j <= 100; j++)
    {
      double p = RandUtil::randUnif(0.0,1.0);
      if( p <= 0.3)
        S.insert(j);
    }

    Sets2.push_back(S);
  }
  vector < int > HS2;
  HSApprox::get_hs_approximation<int>(Sets2, HS2);

  valid = validate_hitting_set(Sets2,HS2);
 
  if(valid)
  {
    cout << "Hitting set for example 2 is valid " << endl;

    for(size_t i = 0; i < HS2.size(); i++)
      cout << HS2[i] << "  ";
    cout << endl;
  }
  else
  {
    cout << "Hitting set for example 2 is not valid " << endl;
    exit(1);
  }
}

bool validate_hs(const vector<Point>& fatP, vector<size_t>& idxs, size_t k,
double epsilon, size_t dim)
{
  assert(idxs.size() <= fatP.size());

  if(fatP.size() == 0)
    return true;

  if(idxs.size() == 0)
    return false;

  assert(idxs[0] >= 0 && idxs[idxs.size() - 1] < fatP.size());

  size_t ndir = RMSUtils::ndir_for_validation(dim);
  vector<Point> random_dirs;
  RMSUtils::get_random_sphere_points(1.0, dim, ndir, random_dirs,true);




  for(size_t i = 0; i < ndir; i++)
  {

    double cor_max = fatP[idxs[0]].dotP(random_dirs[i]);
    for(size_t j = 1; j < idxs.size(); j++)
    {
      double corval = fatP[idxs[j]].dotP(random_dirs[i]);
      if(corval > cor_max)
        cor_max = corval;
    }

    size_t cnt = 0;
    for(size_t j = 0; j < fatP.size(); j++)
    {
      double ptval = fatP[j].dotP(random_dirs[i]);

      if(cor_max < (1 - epsilon) * ptval)
      {
        cnt++;
	if(cnt >= k)
	  return false;
      }
    }
  }

  return true;

}

bool validate_hs2(const vector<Point>& fatP, vector<size_t>& idxs, size_t k, double epsilon, size_t dim, vector<Point>& Sky)
{
  assert(idxs.size() <= fatP.size());

  if(fatP.size() == 0)
    return true;

  if(idxs.size() == 0)
    return false;

  assert(idxs[0] >= 0 && idxs[idxs.size() - 1] < fatP.size());

  size_t ndir = RMSUtils::ndir_for_validation(dim);
  vector<Point> random_dirs;
  RMSUtils::get_random_sphere_points(1.0, dim, ndir, random_dirs,true);




  for(size_t i = 0; i < ndir; i++)
  {

    double cor_max = fatP[idxs[0]].dotP(random_dirs[i]);
    for(size_t j = 1; j < idxs.size(); j++)
    {
      double corval = fatP[idxs[j]].dotP(random_dirs[i]);
      if(corval > cor_max)
        cor_max = corval;
    }

    size_t cnt = 0;
    for(size_t j = 0; j < Sky.size(); j++)
    {
      double ptval = Sky[j].dotP(random_dirs[i]);

      if(cor_max < (1 - epsilon) * ptval)
      {
        cnt++;
	if(cnt >= k)
	  return false;
      }
    }
  }

  return true;

}




double error_hs(const vector<Point>& fatP, vector<size_t>& idxs, size_t k, size_t dim, double &avg, vector<double> & perc){
  
  assert(idxs.size() <= fatP.size());

  if(fatP.size() == 0)
    return true;

  if(idxs.size() == 0)
    return false;

  assert(idxs[0] >= 0 && idxs[idxs.size() - 1] < fatP.size());

  size_t ndir = RMSUtils::ndir_for_validation(dim);
  vector<Point> random_dirs;
  RMSUtils::get_random_sphere_points(1.0, dim, ndir, random_dirs,true);  

  double sum1=0;
  double maxRegret=0;
  double setMax;
  vector<double> Scores;
  for(size_t i = 0; i < ndir; i++)
  {    
      std::priority_queue<double> q;
    double cor_max = fatP[idxs[0]].dotP(random_dirs[i]);
    for(size_t j = 1; j < idxs.size(); j++)
    {
      double corval = fatP[idxs[j]].dotP(random_dirs[i]);
      if(corval > cor_max)
        cor_max = corval;
    }
    
    
    for(size_t j=0;j<fatP.size();j++){
    	q.push(fatP[j].dotP(random_dirs[i]));
    }
	
    for(int j=1;j<=k-1;j++){
    	q.pop();
    }    
    setMax=q.top();
    
    
    
    double eps=0;
    if(cor_max<setMax){
     eps=1-cor_max/setMax;
     if (eps>maxRegret){
      maxRegret=eps;
     }
     Scores.push_back(eps);
     sum1=sum1+eps;
    }else{
      Scores.push_back((double)0);
    }    
    q=priority_queue<double>();
  }
  sort(Scores.begin(), Scores.end());
  int g=(int)(0.8*Scores.size()+0.5);
 // perc80=Scores[g];
  avg=sum1/(double)ndir;
  return(maxRegret);
}


double error_hs2(const vector<Point>& fatP, vector<size_t>& idxs, size_t k, size_t dim, double &avg, vector<double>& perc, const vector<Point>& SkyData ){


  size_t ndir = RMSUtils::ndir_for_validation(dim);
  vector<Point> random_dirs;
  RMSUtils::get_random_sphere_points(1.0, dim, ndir, random_dirs,true);  

  double sum1=0;
  double maxRegret=0;
  vector<double> Scores;
  for(size_t i = 0; i < ndir; i++)
  {    
    double cor_max = fatP[idxs[0]].dotP(random_dirs[i]);
    for(size_t j = 1; j < idxs.size(); j++)
    {
      double corval = fatP[idxs[j]].dotP(random_dirs[i]);
      if(corval > cor_max)
        cor_max = corval;
    }
  
    vector <size_t> tempIDX;
    double setMax;
    for(size_t l=0;l<k;l++){
      size_t idMax;      
      setMax=-1;
      double tempMax;
      for(size_t j = 0; j < SkyData.size(); j++){
	if(std::find(tempIDX.begin(), tempIDX.end(), j) != tempIDX.end()) { //LOOOK AGAIN
	  continue;
	}else{
	  tempMax=SkyData[j].dotP(random_dirs[i]);
	  if (tempMax>setMax){
	    setMax=tempMax;
	    idMax=j;
	  }
	}
      }
      tempIDX.push_back(idMax);
    }    
    double eps=0;
    if(cor_max<setMax){     
     eps=1-cor_max/setMax;
     Scores.push_back(eps);
     sum1=sum1+eps;     
    }else{
      Scores.push_back((double)0);
    }
    if (eps>maxRegret)
      maxRegret=eps;
  }
  sort(Scores.begin(), Scores.end());
  int g=(int)(0.8*Scores.size()+0.5);
  //perc80=Scores[g];
  avg=sum1/(double)ndir;
  return(maxRegret);
}



void hs_by_sampling(vector<Point> fatP, size_t dim, size_t k, double epsilon, size_t ss, vector<size_t>& idxs)
{
  /* saample ss directions on unit sphere */
  vector<Point> randomP;
  RMSUtils::get_random_sphere_points(1.0, dim, ss, randomP, true);

  
  
  /* generate the sets along each such direction */
  vector < set < size_t > > DSets;
  int flag;
  for(size_t i = 0; i < ss; i++)
  {
    set < size_t > DSet;
    vector < size_t > topkI;
    vector < double > topkV;
    /* get the max 10k along the direction randomP[i] */
    size_t k1 = (10 * k <= fatP.size())? 10 * k : fatP.size();
  
    RMSUtils::rank_selection_dotp(fatP, randomP[i], k1 , topkI, topkV);   


    assert(topkV[k1-1] >= 0);
    flag=1;
    for(size_t j=0;j<idxs.size();j++){
      if(fatP[idxs[j]].dotP(randomP[i])>=(1-epsilon)*topkV[k-1]){
	flag=0;
      }
    }
    if(flag==1){
      /* make the set containing the top k indexes */
      for(size_t j = 0; j < k1; j++)
      {
	if(topkV[j] >= (1- epsilon) * topkV[k-1]){
	  DSet.insert(topkI[j]);
	}else{
	  break;
	}
      }
      DSets.push_back(DSet);
    }
  }

  /* use greedy approximation algorithm for hitting set */
  vector<size_t> idxs2;
  HSApprox::get_hs_approximation<size_t>(DSets, idxs2);
  for(size_t j=0;j<idxs2.size();j++){
    idxs.push_back(idxs2[j]);
  }
}


int get_hs(const vector<Point>& fatP, size_t k, double epsilon, vector<Point>& hsP, 
vector<Point> SkyData, int iter, double &time)
{
 
  clockMhs=(float)0.0;
  clock_t clockStart=clock();
  
  if(fatP.size() == 0)
    return 0;

  if(k >= fatP.size())
  {
    hsP.push_back(fatP[0]);
    return 1;
  }

  size_t dim = fatP[0].get_dimension();

  if(fatP.size() < dim)
  {
    hsP = fatP;
    return fatP.size();
  }



  bool found = false;

  
  vector<size_t> idxs;  

  size_t ss = 10;
  idxs.clear();
  double avgEr;
  double perc80;
  while(1)
  {    
      hsP.clear();
    hs_by_sampling(fatP,dim,k,epsilon/5,ss,idxs);    
    for(vector<size_t>::iterator it=idxs.begin();it!=idxs.end();it++)
        hsP.push_back(fatP[*it]);
    double error;
    bool cond;
    if(k>1){
      cond=validate_hs(fatP, idxs, k, epsilon, dim);     
    }else{      
      cond=validate_hs2(fatP, idxs, k, epsilon, dim, SkyData);   
    }

    if(cond){      

	clock_t clockEnd=clock();
	clockMhs=(float)(clockEnd-clockStart) / (CLOCKS_PER_SEC / 1000);
	
        float tTime=clockMhs;
        
        //======================================================================================================
        double error,avgEr;
        vector<double> perc(10);
        if(k>1){	  
	  error=error_hs(fatP, idxs, k, dim, avgEr, perc);          
	}else{
	  error=error_hs2(fatP, idxs, k, dim, avgEr, perc, SkyData);	 
	}        
        //======================================================================================================
        
	// cout<<epsilon<<"\t"<<idxs.size()<<"\t"<<tTime<<"\t"<<error<<"\t"<<avgEr<<"\n"<<std::flush;
    time = tTime;
	return idxs.size();
   
    }
    ss =2*ss;
  }
}

void runHS(vector<Point> dataP, int r, int k,vector<Point> &result, vector<Point> curSky, double &time)
{
    cout << "running HS...." ;
  size_t dim = dataP[0].get_dimension();
  
  vector<Point> hsP;

    double upper = 1.0, lower = 0, epsilon = (upper + lower) /2;
  int core = 0;  
  while(core < r || core > (r+5))
  {
    hsP.clear();
    core = get_hs(dataP, k, epsilon, hsP, curSky,0, time);  
    if(core > (r + 5))
    {
        lower = epsilon;
        epsilon = (epsilon + upper) / 2;
    }
    else 
    {
        upper = epsilon;
        epsilon = (epsilon + lower) / 2;
    }
    if(upper - lower <= 0.001)
        break;
    // cout << core << endl;

  }
  int size = min(r, (int)hsP.size());
  for (size_t i = 0; i < size; i++)
  {
      result.push_back(hsP[i]);
  }
  cout << " " << result.size() << endl;
}


