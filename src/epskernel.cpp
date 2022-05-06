#include"epskernel.h"

ofstream myfile;
float clockM=(float)0.0;
float clockM2=(float)0.0;
double ResR[10];
double ResRe[10];
double ResT[10];


void get_anns(ANN* ann_ds, const vector<Point>& queryP, double
epsilon, vector<size_t>& idxs)
{
  //get the required ANNs
  ann_ds->getANNs(queryP, epsilon, idxs);
}

void get_anns(ANN* ann_ds, const vector<Point>& queryP, double
epsilon, vector<size_t>& idxs, int r)
{
  //get the required ANNs
  ann_ds->getANNks(queryP, epsilon, idxs, r);
}

/** 
    The point set is represented in the vector fatP. The idxs is a *sorted* list
    of indexes into the fatP vector. 
    The function generates a lot of random directions and checks for
    each if the coreset max is at least 1 - epsilon of the max over
    points
    Since the point set is fat the maximum in any direction of the projections
    is +ive.
**/
bool validate_coreset(const vector<Point>& fatP, vector<size_t>& idxs, double epsilon, size_t dim, int k)
{
  //clearly the coreset must be a subset!
  assert(idxs.size() <= fatP.size());
 

  if(fatP.size() == 0)
    return true;

  if(idxs.size() == 0)
    return false;

  assert(idxs[0] >= 0 && idxs[idxs.size() - 1] < fatP.size());

  size_t ndir = RMSUtils::ndir_for_validation(dim);
   
  //cout << "Validate_coreset: ndir = " << ndir << endl; 
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


bool validate_coreset2(const vector<Point>& fatP, vector<size_t>& idxs, double epsilon, size_t dim, int k, const vector<Point>& Skyline)
{
  //clearly the coreset must be a subset!
  assert(idxs.size() <= fatP.size());
 

  if(fatP.size() == 0)
    return true;

  if(idxs.size() == 0)
    return false;

  assert(idxs[0] >= 0 && idxs[idxs.size() - 1] < fatP.size());

  size_t ndir = RMSUtils::ndir_for_validation(dim);
   
  //cout << "Validate_coreset: ndir = " << ndir << endl; 
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
    for(size_t j = 0; j < Skyline.size(); j++)
    {
      double ptval = Skyline[j].dotP(random_dirs[i]);  
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


/**
    Adjoin to the current coreset (represented by idxs) 
    by taking deltam new random points on the outer
    ball (which is known to surround the CH of fatP (see caller fn)) 
    and taking the delta-ANN
    of each such point. The input point set is a fat point set

    The return value is a set of indexes (into fatP) that are unique (i.e.,
    no value is repeated) and sorted.
**/
void coreset_by_sample(ANN* ann_ds, 
                       size_t dim, 
                       double outer_rad, 
                       double delta,
                       size_t deltam,
                       vector<size_t>& idxs, set<size_t> &unique_idxs, int r)
{
  vector<Point> randomP;
  RMSUtils::get_random_sphere_points(outer_rad, dim, deltam, randomP,true);

  //cout << "Returned after getting " << deltam << " random points" << endl;

  //now get the anns of the randomP points on the fat point set
  //The points are stored in vector fatP, the return
  //value will be an array of indexes into the vector fatP of
  //total size randomP
  vector<size_t> new_idxs;
  get_anns(ann_ds, randomP, delta, new_idxs, r);

  //figure out all the distinct indices
  //unordered_set<size_t> unique_idxs;
  
  //for(size_t i = 0; i < idxs.size(); i++)
   // unique_idxs.insert(idxs[i]);
  for(size_t i = 0; i < new_idxs.size(); i++){
    //unique_idxs.insert(new_idxs[i]);
      if(unique_idxs.count(new_idxs[i])<=0){
          idxs.push_back(new_idxs[i]);
          unique_idxs.insert(new_idxs[i]);
      }
  }
  //cout << "# unique indexes =  " << unique_idxs.size() << endl;
  /*
  idxs.clear();

  
  set<size_t>::iterator it;
  //this iteration puts unique indexes in sorted order
  for(it = unique_idxs.begin(); it != unique_idxs.end(); ++it)
    idxs.push_back((*it)); */
}



double error_coreset(const vector<Point>& fatP, vector<size_t>& idxs, size_t k, size_t dim, double &avg, double & perc80){
  
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
  perc80=Scores[g];  
  avg=sum1/(double)ndir;
  return(maxRegret);
}

double error_coreset2(const vector<Point>& fatP, vector<size_t>& idxs, size_t k, size_t dim, double &avg, double & perc80, const vector<Point>& Skyline){
  
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
      for(size_t j = 0; j < Skyline.size(); j++){
	if(std::find(tempIDX.begin(), tempIDX.end(), j) != tempIDX.end()) { //LOOOK AGAIN
	  continue;
	}else{
	  tempMax=Skyline[j].dotP(random_dirs[i]);
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
     if (eps>maxRegret){  
      maxRegret=eps;
     }
     Scores.push_back(eps);
     sum1=sum1+eps;
    }else{
      Scores.push_back((double)0);
    }
  }
  sort(Scores.begin(), Scores.end());
  int g=(int)(0.8*Scores.size()+0.5);
  perc80=Scores[g];  
  avg=sum1/(double)ndir;
  return(maxRegret);
}




int get_coreset(const vector<Point>& fatP, double epsilon, vector<Point>& coresetP, int k, 
const vector<Point> & Skyline, int iter, double &time)
{

  clock_t clockS=clock();
  if(fatP.size() == 0)
    return 0;

  size_t dim = fatP[0].get_dimension();
  //cout << "get_coreset: dim = " << dim << endl;

  if(fatP.size() < dim)
  {
    coresetP = fatP;
    return fatP.size();
  }
  double inner_rad, outer_rad;

  outer_rad=1+sqrt(dim);
  double delta = epsilon / (2 * outer_rad);


  size_t lim = 8 * sizeof(size_t);
  vector<size_t> idxs;
  bool found = false;
  set<size_t> unique_idxs;
  ANN* ann_ds = new ANN();  
  ann_ds->insertPts(fatP);  

  size_t m = 10;

  while(1)
  {      
      coresetP.clear();
    coreset_by_sample(ann_ds, dim, outer_rad, delta, m, idxs, unique_idxs, k);
    for(set<size_t>::iterator it=unique_idxs.begin();it!=unique_idxs.end();it++)
        coresetP.push_back(fatP[*it]);
    bool b;
    if(k>1){
      b=validate_coreset(fatP, idxs, epsilon, dim,k);   
    }else{
      b=validate_coreset2(fatP, idxs, epsilon, dim,k, Skyline);  
    }
    if(b)
    {
      found = true;            
      goto done;
    }    
    m *= 2;
  } 


done:

  delete ann_ds;
  double avgEr;
  double perc80;
  if(found)
  {
  

    double Merror;
    clock_t clockE=clock();
    float tTime=(float)(clockE-clockS) / (CLOCKS_PER_SEC / 1000);

    if(k>1){
      Merror=error_coreset(fatP, idxs, k, dim, avgEr, perc80);  
    }else{
      Merror=error_coreset2(fatP, idxs, k, dim, avgEr, perc80, Skyline);   
    }


    myfile<<"\n";
    myfile<<epsilon<<"\t"<<idxs.size()<<"\t"<<Merror<<"\t"<<avgEr<<"\t"<<perc80<<"\t"<<tTime<<"\n";
    // cout<<epsilon<<"\t"<<idxs.size()<<"\t"<<Merror<<"\t"<<avgEr<<"\t"<<perc80<<"\t"<<tTime<<"\n"<<std::flush;    
    time = tTime;
    ResR[iter]=ResR[iter]+idxs.size();
    ResRe[iter]=ResRe[iter]+Merror;
    ResT[iter]=ResT[iter]+tTime;
    return idxs.size();
  }
  else
    return -1;
}



void runEpsKernel(vector<Point> dataP, int r, int k,vector<Point> &result, vector<Point> curSky, double &time)
{
    cout << "running Eps-Kernel...." ;
  size_t dim = dataP[0].get_dimension();
  
  vector<Point> coresetP;

    double upper = 1.0, lower = 0, epsilon = (upper + lower) /2;
  int core = 0;  
  while(core < r || core > (r+5))
  {
    coresetP.clear();
    core = get_coreset(dataP, epsilon, coresetP, k, curSky,0, time);  
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
  int size = min(r, (int)coresetP.size());
  for (size_t i = 0; i < size; i++)
  {
      result.push_back(coresetP[i]);
  }
  cout << " " << result.size() << endl;
}
