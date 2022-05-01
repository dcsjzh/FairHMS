#include "BiGreedy.h"

/* compute the increased benefits and sort */
class Increment
{
public:
      Increment();
      Increment(int i, double x)
      {
            this->value = x;
            this->index = i;
      }
      bool operator<(const Increment &right) const { return value < right.value; }
      // bool operator<(const Increment &right) const { return abs(value - right.value) < 1e-7 ? index > right.index : value < right.value; }
      int get_index() const { return index; }
      double get_value() const { return value; }

private:
      int index;
      double value;
};


/* store the increased benefits with a priority queue */
inline void updatePriorityQueue(vector<Point> &dataP, vector<int> &candidate, vector<int> &result, vector<Point> &allUtility,
                                int m, double epsilon, vector<pair<int, double>> &maxD, vector<pair<int, double>> &maxS, double lastMaxDelta,
                                priority_queue<Increment> &prque)
{
      while (!prque.empty())
      {
            int pid = prque.top().get_index();
            prque.pop();
            if (find(candidate.begin(), candidate.end(), pid) == candidate.end())
            {
                  continue;
            }
            Point p = dataP[pid];
            double tmpDelta = 0.0;
            for (int i = 0; i < m; ++i)
            {
                  double newValue = p.dotP(allUtility[i]);
                  if (newValue > maxS[i].second)
                  {
                        tmpDelta += min(newValue / maxD[i].second, 1.0 - epsilon) - min(maxS[i].second / maxD[i].second, 1.0 - epsilon);
                  }
            }
            prque.push(Increment(p.getId(), tmpDelta / m));
            if (abs(tmpDelta - lastMaxDelta * m) <= 1e-9)
            {
                  //Remove points not in candidate
                  break;
            }
            lastMaxDelta = prque.top().get_value();
      }
}

/* construct the candidate set */
vector<int> constructCandidate(unordered_map<int, vector<Point>> &categorizedDataP, 
                                      unordered_map<int, fairinCate> &partitionMatroid,
                                      vector<int> &result,
                                      unordered_map<int, int> fairofR,
                                      int k)
{
      vector<int> candidate;
      int remained = k - result.size(); 
      for (auto it = partitionMatroid.begin(); it != partitionMatroid.end(); it++)
      {
            int cateId = it->first;
            if (it->second.lc > fairofR[cateId])
            {
                  for (Point p : categorizedDataP[cateId])
                  {
                        if (find(result.begin(), result.end(), p.getId()) == result.end())
                        {
                              candidate.push_back(p.getId());
                        }
                  }
                  remained -= it->second.lc;
            }
      }
      for (auto it = partitionMatroid.begin(); it != partitionMatroid.end(); it++)
      {
            int cateId = it->first;
            if (it->second.lc <= fairofR[cateId] && it->second.uc > fairofR[cateId] && remained > 0)
            {
                  for (Point p : categorizedDataP[cateId])
                  {
                        if (find(result.begin(), result.end(), p.getId()) == result.end())
                        {
                              candidate.push_back(p.getId());
                        }
                  }
            }
      }
      return candidate;
}

/* the Multi-round Greedy algorithm */
void multiroundGreedy(unordered_map<int, vector<Point>> &categorizedDataP, 
                      unordered_map<int, fairinCate> &partitionMatroid,   
                      vector<Point> &dataP,                               
                      int cateId,                                         
                      double epsilon,                                     
                      vector<Point> &allUtility,   
                      vector<int> &result,
                      vector<pair<int, double>> &maxD, vector<pair<int, double>> &maxS,
                      double tau,
                      unordered_map<int, int> &fairofR,
                      int k, int m, int gamma,
                      double &totalDelta)
{
      result.clear();
      for (auto it = fairofR.begin(); it != fairofR.end(); it++)
      {
            it->second = 0;
      }
      for (int i = 0; i < m; ++i)
      {
            maxS[i].first = -1;
            maxS[i].second = 0.0;
      }
      vector<int> tmpResult;
      vector<pair<int, double>> maxTmpS = maxS, maxTmpD = maxD;
      for (int count = 0; count < gamma; ++count)
      {
            tmpResult.clear();
            for (int i = 0; i < m; ++i)
            {
                  maxTmpS[i].first = -1;
                  maxTmpS[i].second = 0.0;
            }
            vector<int> candidate = constructCandidate(categorizedDataP, partitionMatroid, result, fairofR, k);
            
            priority_queue<Increment> prque;
            for (int pid : candidate)
            {
                  Point p = dataP[pid];
                  double tmpDelta = 0.0;
                  for (int i = 0; i < m; ++i)
                  {
                        double newValue = p.dotP(allUtility[i]);
                        if (newValue > maxTmpS[i].second)
                        {
                              tmpDelta += min(newValue / maxTmpD[i].second, tau) - min(maxTmpS[i].second / maxTmpD[i].second, tau);
                        }
                  }
                  prque.push(Increment(p.getId(), tmpDelta / m));
            }
            
            while (candidate.size() > 0 && result.size() < k && tmpResult.size() < k)
            {
                  int id = -1;
                  double maxDelta = 0.0;
                  do
                  {
                        id = prque.top().get_index();
                        maxDelta = prque.top().get_value();
                        prque.pop();
                  }
                  while (find(candidate.begin(), candidate.end(), id) == candidate.end());

                  if (maxDelta - 0.0 < 1e-9)
                  {
                        break;
                  }
                  Point p = dataP[id];
                  //recalculate the increased benefits with maxD
                  double tmpDelta = 0.0;
                  for (int i = 0; i < m; ++i)
                  {
                        double newValue = p.dotP(allUtility[i]);
                        if (newValue > maxS[i].second)
                        {
                              tmpDelta += min(newValue / maxD[i].second, tau) - min(maxS[i].second / maxD[i].second, tau);
                        }
                  }
                  totalDelta += tmpDelta / m;
                  for (int i = 0; i < m; ++i)
                  {
                        maxTmpS[i].first = p.getId();
                        maxTmpS[i].second = max(maxTmpS[i].second, p.dotP(allUtility[i]));
                        maxS[i].first = p.getId();
                        maxS[i].second = max(maxS[i].second, p.dotP(allUtility[i]));
                  }
                  
                  tmpResult.push_back(id);
                  result.push_back(id);
                  ++fairofR[p.get_category(cateId)];
                  candidate = constructCandidate(categorizedDataP, partitionMatroid, result, fairofR, k);
                  updatePriorityQueue(dataP, candidate, tmpResult, allUtility, m, 1 - tau, maxTmpD, maxTmpS, maxDelta, prque);
            }
            if (tmpResult.size() > 0)
            {
                  for (int i = 0; i < m; ++i)
                  {
                        double maxUofD = 0.0;
                        int maxIdx = -1;
                        for (int j = 0; j < dataP.size(); ++j)
                        {
                              if (find(tmpResult.begin(), tmpResult.end(), j) != tmpResult.end())
                              {
                                    continue;
                              }
                              if (maxUofD < dataP[j].dotP(allUtility[i]))
                              {
                                    maxUofD = dataP[j].dotP(allUtility[i]);
                                    maxIdx = j;
                              }
                        }
                        maxTmpD[i].first = maxIdx;
                        maxTmpD[i].second = maxUofD;
                  }
            }
            if (totalDelta >= (1 - epsilon / (2 * m)) * tau)
            {
                  break;
            }
      }
      if (totalDelta < (1 - epsilon / (2 * m)) * tau)
      {
            result.clear();
            totalDelta = 0;
            for (auto it = fairofR.begin(); it != fairofR.end(); it++)
            {
                  it->second = 0;
            }
      }
}

/* use binary search to find the maximum \tau */
void binSearchforResult(unordered_map<int, vector<Point>> &categorizedDataP,
                        unordered_map<int, fairinCate> &partitionMatroid,  
                        vector<Point> &dataP,                              
                        int cateId,                                        
                        double epsilon,                                    
                        vector<Point> &allUtility,                        
                        vector<int> &result,
                        vector<pair<int, double>> &maxD, vector<pair<int, double>> &maxS,
                        double tauHigh, double tauLow, double &tau, 
                        unordered_map<int, int> &fairofR,
                        int k, int m, int gamma)
{
      double curDelta = 0.0;
      while (tauHigh - tauLow - epsilon >= 1e-9)
      {
            tau = (tauHigh + tauLow) / 2;
            double totalDelta = 0.0;
            vector<int> tmpResult;
            unordered_map<int, int> tmpFairofR;
            for (auto it = fairofR.begin(); it != fairofR.end(); ++it)
            {
                  tmpFairofR[it->first] = 0;
            }
            multiroundGreedy(categorizedDataP, partitionMatroid, dataP, cateId, epsilon, allUtility, tmpResult,
                             maxD, maxS, tau, tmpFairofR, k, m, gamma, totalDelta);
            if (tmpResult.size() != 0 && totalDelta > curDelta)
            {
                  result = tmpResult;
                  fairofR = tmpFairofR;
                  curDelta = totalDelta;
                  tauLow = tau;
            }
            else
            {
                  tauHigh = tau;
            }
      }
}

void runBiGreedyAlg(unordered_map<int, vector<Point>> &categorizedDataP, 
                      unordered_map<int, fairinCate> &partitionMatroid,   
                      vector<Point> &dataP,                               
                      int cateId,                                         
                      double epsilon,                                     
                      vector<Point> &allUtility,                          
                      vector<int> &result,
                      int k, int maxM, double &time)
{
      int numofCategory = categorizedDataP.size();
      unordered_map<int, int> fairofR, tmpFairofR;
      
      for (auto it = partitionMatroid.begin(); it != partitionMatroid.end(); ++it)
      {
            fairofR[it->first] = 0;
      }
      int dim = dataP[0].get_dimension();
      int m; 
      m = k * dim * 5; 
      m = min(m, maxM);
      
      vector<pair<int, double>> maxD(m, make_pair(-1, 0.0)), maxS(m, make_pair(-1, 0.0)); //int-id double-utility
      for (int i = 0; i < m; ++i)
      {
            double maxUofD = 0.0;
            int maxIdx = -1;
            for (Point p : dataP)
            {
                  int pid = p.getId();
                  if (maxUofD < p.dotP(allUtility[i]))
                  {
                        maxUofD = p.dotP(allUtility[i]);
                        maxIdx = pid;
                  }
            }
            maxD[i].first = maxIdx;
            maxD[i].second = maxUofD;
      }
      
      double tauHigh = 1.0, tauLow = 0.0, tau;

      clock_t start = clock();
      double up2now = 0.0; 
      int gamma = log(2 * m / epsilon) / log(2);
      binSearchforResult(categorizedDataP, partitionMatroid, dataP, cateId, epsilon, allUtility,
                         result, maxD, maxS, tauHigh, tauLow, tau, fairofR, k, m, gamma);
                         
      clock_t end = clock();
      time = (double)(end - start) / (CLOCKS_PER_SEC / 1000);
      vector<int> candidate = constructCandidate(categorizedDataP, partitionMatroid, result, fairofR, k);
      
      while (candidate.size() > 0 && result.size() < k)
      {
            int id = candidate[0];
            Point p = dataP[id];
            result.push_back(id);
            for (int i = 0; i < m; ++i)
            {
                  maxS[i].first = p.getId();
                  maxS[i].second = max(maxS[i].second, p.dotP(allUtility[i]));
            }
            ++fairofR[p.get_category(cateId)];
            candidate = constructCandidate(categorizedDataP, partitionMatroid, result, fairofR, k);
      }
}

void runBiGreedyPlusAlg(unordered_map<int, vector<Point>> &categorizedDataP, 
                          unordered_map<int, fairinCate> &partitionMatroid,   
                          vector<Point> &dataP,                               
                          int cateId,                                         
                          double lambda,                                      
                          double epsilon,                                     
                          vector<Point> &allUtility,                          
                          vector<int> &result,
                          int k, int maxM, double &time)
{
      int numofCategory = categorizedDataP.size();
      unordered_map<int, int> fairofR; 
      vector<int> tmpResult;
      unordered_map<int, int> fairofTmpR;
      
	  for (auto it = partitionMatroid.begin(); it != partitionMatroid.end(); ++it)
      {
            fairofR[it->first] = 0;
            fairofTmpR[it->first] = 0;
      }
      int dim = dataP[0].get_dimension();
      int m; 
      m = dim * k * 1; 
      m = min(m, maxM);
      
      vector<pair<int, double>> maxD(m, make_pair(-1, 0.0)), maxS(m, make_pair(-1, 0.0)); //int-id double-utility
      for (int i = 0; i < m; ++i)
      {
            double maxUofD = 0.0;
            int maxIdx = -1;
            for (Point p : dataP)
            {
                  int pid = p.getId();
                  if (maxUofD < p.dotP(allUtility[i]))
                  {
                        maxUofD = p.dotP(allUtility[i]);
                        maxIdx = pid;
                  }
            }
            maxD[i].first = maxIdx;
            maxD[i].second = maxUofD;
      }
      
      double up2now = 0.0; 
      double tauHigh = 1.0, tauLow = 0.0, tau, preTau;
      int gamma = log(2 * m / epsilon) / log(2);
      clock_t start = clock();
      binSearchforResult(categorizedDataP, partitionMatroid, dataP, cateId, epsilon, allUtility,
                         result, maxD, maxS, tauHigh, tauLow, tau, fairofR, k, m, gamma);
      while (true)
      {
            m *= 2;
            m  = min(m, maxM);
            for (int i = m / 2; i < m; ++i)
            {
                  double maxUofD = 0.0;
                  int maxIdx = -1;
                  for (Point p : dataP)
                  {
                        int pid = p.getId();
                        if (maxUofD < p.dotP(allUtility[i]))
                        {
                              maxUofD = p.dotP(allUtility[i]);
                              maxIdx = pid;
                        }
                  }
                  maxD.push_back(make_pair(maxIdx, maxUofD));
                  maxS.push_back(make_pair(-1, 0.0));
            }
            vector<pair<int, double>> maxTmpS(m, make_pair(-1, 0.0));
            double totalDelta;
            multiroundGreedy(categorizedDataP, partitionMatroid, dataP, cateId, epsilon, allUtility, tmpResult,
                             maxD, maxTmpS, tau - lambda, fairofTmpR, k, m, gamma, totalDelta);
            if (tmpResult.size() != 0)
            {
                  binSearchforResult(categorizedDataP, partitionMatroid, dataP, cateId, epsilon, allUtility,
                                     result, maxD, maxS, tau, tau - lambda, tau, fairofR, k, m, gamma);
                  break;
            }
            else
            {
                  binSearchforResult(categorizedDataP, partitionMatroid, dataP, cateId, epsilon, allUtility,
                                     result, maxD, maxS, tau - lambda, 0, tau, fairofR, k, m, gamma);
            }
      }
      clock_t end = clock();
      time = (double)(end - start) / (CLOCKS_PER_SEC / 1000);
      vector<int> candidate = constructCandidate(categorizedDataP, partitionMatroid, result, fairofR, k);
      while (candidate.size() > 0 && result.size() < k)
      {
            int id = candidate[0];
            Point p = dataP[id];
            result.push_back(id);
            for (int i = 0; i < m; ++i)
            {
                  maxS[i].first = p.getId();
                  maxS[i].second = max(maxS[i].second, p.dotP(allUtility[i]));
            }
            ++fairofR[p.get_category(cateId)];
            candidate = constructCandidate(categorizedDataP, partitionMatroid, result, fairofR, k);
      }
}