#include "greedy.h"

size_t returnMax1Coordinate(vector<Point> &fatP)
{
  double maxC;
  size_t I;
  maxC = fatP[0].get_coordinate(0);
  I = 0;
  for (int i = 1; i < fatP.size(); i++)
  {
    if (fatP[i].get_coordinate(0) > maxC)
    {
      maxC = fatP[i].get_coordinate(0);
      I = i;
    }
  }
  return I;
}

void runGreedy(vector<Point> dataP, int r, int k, vector<Point> &R, vector<Point> curSky, double &time)
{
  cout << "Running Greedy..." << endl;
  clock_t tstart = clock();

  std::srand(std::time(0));

  GRBEnv env = GRBEnv();

  size_t dim = dataP[0].get_dimension();

  vector<Point> fatP = dataP;
  assert(r >= 1);

  size_t ndir = RMSUtils::ndir_for_validation(dim);

  vector<size_t> idxs;

  size_t maxIndex = returnMax1Coordinate(fatP);
  idxs.push_back(maxIndex);
  R.push_back(fatP[maxIndex]);

  double maxRegret = 0.0;
  double MaxRegret_p;
  int optimstatus;
  size_t id;
  double STime = (double)0;
  clock_t tend = clock();
  double tTime = (double)(tend - tstart) / (CLOCKS_PER_SEC / 1000);
  STime += tTime;
  while (idxs.size() < r)
  {
    clock_t tStart_a = clock();
    maxRegret = 0;
    //cout<<"================================================= "<<idxs.size()<<"\n";
    for (size_t i = 0; i < fatP.size(); i++)
    {
      if (find(idxs.begin(), idxs.end(), i) != idxs.end())
      {
        continue;
      }
      else
      {
        GRBModel model = GRBModel(env);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        GRBVar X[dim + 1];
        for (size_t j = 0; j < dim + 1; j++)
        {
          X[j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }
        model.update();
        GRBLinExpr exprObj = X[dim];
        model.setObjective(exprObj, GRB_MAXIMIZE);

        //Equalities
        GRBLinExpr exprEq = GRBLinExpr();
        exprEq = fatP[i].get_coordinate(0) * X[0];
        for (size_t j = 1; j < dim; j++)
        {
          exprEq = exprEq + fatP[i].get_coordinate(j) * X[j];
          //exprEq.addTerms(fatP[i].get_coordinate(j),X[j],0);
        }
        model.addConstr(exprEq == 1.0);

        //Inequalities
        for (size_t j = 0; j < R.size(); j++)
        {
          GRBLinExpr exprIneq = GRBLinExpr();
          exprIneq = -1.0 * X[dim];
          for (size_t z = 0; z < dim; z++)
          {
            exprIneq = exprIneq + (fatP[i].get_coordinate(z) - R[j].get_coordinate(z)) * X[z];
          }
          model.addConstr(exprIneq >= 0);
        }
        model.update();
        model.optimize();
        optimstatus = model.get(GRB_IntAttr_Status);
        if (optimstatus == GRB_OPTIMAL)
        {
          MaxRegret_p = X[dim].get(GRB_DoubleAttr_X);
          if (MaxRegret_p > maxRegret)
          {
            maxRegret = MaxRegret_p;
            id = i;
          }
        }
        // model.reset();
      } //end else
    }

    clock_t tStart_b = clock();
    double tTime = (double)(tStart_b - tStart_a) / (CLOCKS_PER_SEC / 1000);
    STime = STime + tTime;
    // double MaxR;
    // double AvgR;
    // double perc80;
    // RMSUtils::Max_Avg_Regret(1.0, dim, ndir, fatP, R, MaxR, AvgR, perc80, 1);

    // // cout<<idxs.size()<<"\t"<<maxRegret<<"\t"<<MaxR<<"\t"<<AvgR<<"\t"<<perc80<<"\t"<<tTime<<"\t"<<STime<<"\n"<<std::flush;
    // if (maxRegret == 0.0)
    //   break;

    idxs.push_back(id);
    R.push_back(fatP[id]);
  }

  time = STime;
  clock_t end = clock();
  time = (double)(end - tstart) / (CLOCKS_PER_SEC / 1000);
}
