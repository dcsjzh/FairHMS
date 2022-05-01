#include "greedyK.h"

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

void runGreedyK(vector<Point> dataP, int r, int k, vector<Point> &R, vector<Point> curSky, double &time)
{
  cout << "Running GreedyK...";
  GRBEnv env = GRBEnv();

  size_t dim = dataP[0].get_dimension();
  vector<Point> fatP = curSky;

  assert(r >= 1);

  if (k < 2)
  {
    cout << "k should be greater than 2\n";
    exit(0);
  }

  int T;
  T = 55;

  size_t ndir = RMSUtils::ndir_for_validation(dim); 

  vector<size_t> idxs;

  clock_t tStart1 = clock();
  size_t maxIndex = returnMax1Coordinate(fatP);
  idxs.push_back(maxIndex);    
  R.push_back(fatP[maxIndex]); 
  clock_t tStart2 = clock();

  double MaxR;
  double AvgR;
  double perc80;

  RMSUtils::fast_Max_Avg_Regret(1.0, dim, ndir, dataP, R, MaxR, AvgR, perc80, k); 
                                                                                  //   cout<<"1\t"<<MaxR<<"\t"<<AvgR<<"\t"<<perc80<<"\t"<<((double)(tStart2 - tStart1)/CLOCKS_PER_SEC)<<"----Choose id: "<<maxIndex<<std::flush;

  cout << "\n"
       << std::flush;

  double maxRegret = -1.0;
  double MaxRegret_p;
  int optimstatus;
  double err = 0.000001;
  size_t id;
  int flag1;
  double STime = (double)0;
  while (idxs.size() < r)
  {

    clock_t tStart_a = clock(); 
    maxRegret = -1.0;

    flag1 = 0;
    vector<double> W;
    double sumH = 0.0;
    int countH = 0;
    for (size_t i = 0; i < fatP.size(); i++)
    {
      if (find(idxs.begin(), idxs.end(), i) != idxs.end())
      {
        continue; 
      }
      else
      {
        countH = countH + 1;
        int h;
        for (h = 1; h <= T; h++)
        {
          vector<size_t> Partition[k - 1];
          vector<size_t> Temp;
          for (size_t z = 0; z < fatP.size(); z++)
          {
            if (find(idxs.begin(), idxs.end(), z) != idxs.end())
            {
              continue;
            }
            else
            {
              if (z != i)
              {
                Temp.push_back(z); 
              }
            }
          }
          std::random_shuffle(Temp.begin(), Temp.end()); 

          int group = (int)Temp.size() / (int)(k - 1);
          int z2 = 0;
          for (int z1 = 0; z1 < k - 1; z1++)
          {
            while (z2 < (z1 + 1) * (int)group)
            {
              Partition[z1].push_back(Temp[z2]); 
              z2 = z2 + 1;
            }
          }
          for (int z1 = z2; z1 < Temp.size(); z1++)
          {
            int l = std::rand() % (k - 1);
            Partition[l].push_back(Temp[z1]); 
          }

          GRBModel model = GRBModel(env);
          model.getEnv().set(GRB_IntParam_OutputFlag, 0);
          GRBVar X[dim + k];
          for (size_t j = 0; j < dim + 1; j++)
          {
            int a = j;
            string sname = "X";
            stringstream ss;
            ss << a;
            ss >> sname;
            X[j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, sname); 
          }
          for (size_t j = dim + 1; j < dim + k; j++)
          {
            int a = j;
            string sname = "X";
            stringstream ss;
            ss << a;
            ss >> sname;
            X[j] = model.addVar(-err, GRB_INFINITY, 0.0, GRB_CONTINUOUS, sname);
          }
          model.update();
          GRBLinExpr exprObj = X[dim]; 
          for (size_t j = dim + 1; j < dim + k; j++)
          {
            exprObj = exprObj - err * X[j]; 
  		}
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
          for (size_t j = 0; j < k - 1; j++)
          {
            for (size_t z = 0; z < Partition[j].size(); z++)
            {
              GRBLinExpr exprIneqK = GRBLinExpr();
              exprIneqK = -1.0 * X[dim + j + 1];
              for (size_t zd = 0; zd < dim; zd++)
              {
                exprIneqK = exprIneqK + (fatP[(Partition[j])[z]].get_coordinate(zd) - fatP[i].get_coordinate(zd)) * X[zd];
              }
              model.addConstr(exprIneqK <= 0);
            }
          }
          model.update();

          model.optimize();
         
          optimstatus = model.get(GRB_IntAttr_Status);
          if (optimstatus == GRB_OPTIMAL)
          {
            int flagX = 1;
            for (size_t j = dim + 1; j < dim + k; j++)
            { 
              if (X[j].get(GRB_DoubleAttr_X) <= 0)
              {
                flagX = 0;
                break;
              }
            }

            if (flagX == 1)
            {
              MaxRegret_p = X[dim].get(GRB_DoubleAttr_X); 
              if (MaxRegret_p > maxRegret)
              {
                maxRegret = MaxRegret_p;
                id = i;
                flag1 = 1;
                W.clear();
                for (size_t j = 0; j < dim; j++)
                {
                  W.push_back(X[j].get(GRB_DoubleAttr_X)); 
                }
              }
              break;
            }
          }
          // model.reset();
        } 
        sumH = sumH + (double)h;
      } 
    }  
    if (flag1 == 1)
    {
      vector<size_t> Vtemp;
      Vtemp.push_back(id); 
      for (size_t i = 0; i < fatP.size(); i++)
      {
        if (find(idxs.begin(), idxs.end(), i) != idxs.end())
        {
          continue;
        }
        else
        {
          if (i != id)
          {
            Point Pw = Point(dim, W); 
            if (fatP[i].dotP(Pw) >= fatP[id].dotP(Pw))
            {
              Vtemp.push_back(i);
            }
          }
        }
      }
      id = -1;
      double sumP = 0.0;
      double tempSum;
      for (size_t j = 0; j < Vtemp.size(); j++)
      {
        tempSum = 0.0;
        for (size_t z = 0; z < dim; z++)
        {
          tempSum = tempSum + fatP[Vtemp[j]].get_coordinate(z); 
        }
        if (tempSum > sumP)
        {
          sumP = tempSum;
          id = Vtemp[j];
        }
      }
      idxs.push_back(id);
      R.push_back(fatP[id]);
    }
    else
    { 
      //   cout<<"Point from here (v2) "<<maxRegret<<"\n"<<std::flush;
      id = -1;
      double sumP = 0.0;
      double tempSum;
      for (size_t i = 0; i < fatP.size(); i++)
      {
        if (find(idxs.begin(), idxs.end(), i) != idxs.end())
        {
          continue;
        }
        else
        {
          tempSum = 0.0;
          for (size_t z = 0; z < dim; z++)
          {
            tempSum = tempSum + fatP[i].get_coordinate(z);
          }
          if (tempSum > sumP)
          {
            sumP = tempSum;
            id = i;
          }
        }
      }
      idxs.push_back(id);
      R.push_back(fatP[id]);
    }
    clock_t tStart_b = clock(); 
    double tTime = ((double)(tStart_b - tStart_a) / (CLOCKS_PER_SEC / 1000));
    STime = STime + tTime;
    RMSUtils::fast_Max_Avg_Regret(1.0, dim, ndir, dataP, R, MaxR, AvgR, perc80, k);
    // cout<<idxs.size()<<"\t"<<MaxR<<"\t"<<AvgR<<"\t"<<perc80<<"\t"<<tTime<<"\t"<<" || "<<std::flush;
    // cout<<"---Choose id:"<<id<<"----- TimeSoFar: "<<STime<<"\n"<<std::flush;
  }

  
  //   for(size_t i=0;i<idxs.size();i++){
  //     cout<<"\n"<<std::flush;
  //     for(size_t j=0;j<dim;j++){
  //       cout<<R[i].get_coordinate(j)<<" "<<std::flush;
  //     }
  //   }

  time = STime;
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
