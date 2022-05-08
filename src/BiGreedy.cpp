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
inline void updatePriorityQueue(vector<Point> &dataP, vector<int> &candidate, vector<int> &result, vector<Point> &utiFunClass,
                                int m, double epsilon, vector<pair<int, double>> &maxUtilityofD, vector<pair<int, double>> &maxUtilityofS, double lastMaxDelta,
                                priority_queue<Increment> &priQueofInc)
{
      while (!priQueofInc.empty())
      {
            int pIdx = priQueofInc.top().get_index();
            priQueofInc.pop();
            if (find(candidate.begin(), candidate.end(), pIdx) == candidate.end())
            {
                  continue;
            }
            Point p = dataP[pIdx];
            double tmpDelta = 0.0;
            for (int i = 0; i < m; ++i)
            {
                  double tmpUtility = p.dotP(utiFunClass[i]);
                  if (tmpUtility > maxUtilityofS[i].second)
                  {
                        tmpDelta += min(tmpUtility / maxUtilityofD[i].second, 1.0 - epsilon) - min(maxUtilityofS[i].second / maxUtilityofD[i].second, 1.0 - epsilon);
                  }
            }
            priQueofInc.push(Increment(p.getId(), tmpDelta / m));
            if (abs(tmpDelta - lastMaxDelta * m) <= 1e-9)
            {
                  break;
            }
            lastMaxDelta = priQueofInc.top().get_value();
      }
}

/* construct the candidate set */
vector<int> constructCandidate(unordered_map<int, vector<Point>> &groupedDataP, 
                                      unordered_map<int, fairofGroup> &fairnessConstraint,
                                      vector<int> &result,
                                      unordered_map<int, int> fairofResult,
                                      int k)
{
      vector<int> candidate;
      int remaining = k - result.size(); 
      for (auto it = fairnessConstraint.begin(); it != fairnessConstraint.end(); it++)
      {
            int groupID = it->first;
            if (it->second.lc > fairofResult[groupID])
            {
                  for (Point p : groupedDataP[groupID])
                  {
                        if (find(result.begin(), result.end(), p.getId()) == result.end())
                        {
                              candidate.push_back(p.getId());
                        }
                  }
                  remaining -= it->second.lc;
            }
      }
      for (auto it = fairnessConstraint.begin(); it != fairnessConstraint.end(); it++)
      {
            int groupID = it->first;
            if (it->second.lc <= fairofResult[groupID] && it->second.uc > fairofResult[groupID] && remaining > 0)
            {
                  for (Point p : groupedDataP[groupID])
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
void multiroundGreedy(unordered_map<int, vector<Point>> &groupedDataP, 
                      unordered_map<int, fairofGroup> &fairnessConstraint,   
                      vector<Point> &dataP,                               
                      int groupID,                                         
                      double epsilon,                                     
                      vector<Point> &utiFunClass,   
                      vector<int> &result,
                      vector<pair<int, double>> &maxUtilityofD, vector<pair<int, double>> &maxUtilityofS,
                      double tau,
                      unordered_map<int, int> &fairofResult,
                      int k, int m, int gamma,
                      double &sumDelta)
{
      result.clear();
      for (auto it = fairofResult.begin(); it != fairofResult.end(); it++)
      {
            it->second = 0;
      }
      for (int i = 0; i < m; ++i)
      {
            maxUtilityofS[i].first = -1;
            maxUtilityofS[i].second = 0.0;
      }
      vector<int> tmpResult;
      vector<pair<int, double>> tmpMaxS = maxUtilityofS, tmpMaxD = maxUtilityofD;
      for (int count = 0; count < gamma; ++count)
      {
            tmpResult.clear();
            for (int i = 0; i < m; ++i)
            {
                  tmpMaxS[i].first = -1;
                  tmpMaxS[i].second = 0.0;
            }
            vector<int> candidate = constructCandidate(groupedDataP, fairnessConstraint, result, fairofResult, k);
            
            priority_queue<Increment> priQueofInc;
            for (int pIdx : candidate)
            {
                  Point p = dataP[pIdx];
                  double tmpDelta = 0.0;
                  for (int i = 0; i < m; ++i)
                  {
                        double tmpUtility = p.dotP(utiFunClass[i]);
                        if (tmpUtility > tmpMaxS[i].second)
                        {
                              tmpDelta += min(tmpUtility / tmpMaxD[i].second, tau) - min(tmpMaxS[i].second / tmpMaxD[i].second, tau);
                        }
                  }
                  priQueofInc.push(Increment(p.getId(), tmpDelta / m));
            }
            
            while (candidate.size() > 0 && result.size() < k && tmpResult.size() < k)
            {
                  int pIdx = -1;
                  double maxDelta = 0.0;
                  do
                  {
                        pIdx = priQueofInc.top().get_index();
                        maxDelta = priQueofInc.top().get_value();
                        priQueofInc.pop();
                  }
                  while (find(candidate.begin(), candidate.end(), pIdx) == candidate.end());

                  if (maxDelta - 0.0 < 1e-9)
                  {
                        break;
                  }
                  Point p = dataP[pIdx];
                  //recalculate the increased benefits
                  double tmpDelta = 0.0;
                  for (int i = 0; i < m; ++i)
                  {
                        double tmpUtility = p.dotP(utiFunClass[i]);
                        if (tmpUtility > maxUtilityofS[i].second)
                        {
                              tmpDelta += min(tmpUtility / maxUtilityofD[i].second, tau) - min(maxUtilityofS[i].second / maxUtilityofD[i].second, tau);
                        }
                  }
                  sumDelta += tmpDelta / m;
                  for (int i = 0; i < m; ++i)
                  {
                        tmpMaxS[i].first = p.getId();
                        tmpMaxS[i].second = max(tmpMaxS[i].second, p.dotP(utiFunClass[i]));
                        maxUtilityofS[i].first = p.getId();
                        maxUtilityofS[i].second = max(maxUtilityofS[i].second, p.dotP(utiFunClass[i]));
                  }
                  
                  tmpResult.push_back(pIdx);
                  result.push_back(pIdx);
                  ++fairofResult[p.get_category(groupID)];
                  candidate = constructCandidate(groupedDataP, fairnessConstraint, result, fairofResult, k);
                  updatePriorityQueue(dataP, candidate, tmpResult, utiFunClass, m, 1 - tau, tmpMaxD, tmpMaxS, maxDelta, priQueofInc);
            }
            if (tmpResult.size() > 0)
            {
                  for (int i = 0; i < m; ++i)
                  {
                        double maxUtility = 0.0;
                        int maxIdx = -1;
                        for (int j = 0; j < dataP.size(); ++j)
                        {
                              if (find(tmpResult.begin(), tmpResult.end(), j) != tmpResult.end())
                              {
                                    continue;
                              }
                              if (maxUtility < dataP[j].dotP(utiFunClass[i]))
                              {
                                    maxUtility = dataP[j].dotP(utiFunClass[i]);
                                    maxIdx = j;
                              }
                        }
                        tmpMaxD[i].first = maxIdx;
                        tmpMaxD[i].second = maxUtility;
                  }
            }
            if (sumDelta >= (1 - epsilon / (2 * m)) * tau)
            {
                  break;
            }
      }
      if (sumDelta < (1 - epsilon / (2 * m)) * tau)
      {
            result.clear();
            sumDelta = 0;
            for (auto it = fairofResult.begin(); it != fairofResult.end(); it++)
            {
                  it->second = 0;
            }
      }
}

/* use binary search to find the maximum \tau */
void biSearchResult(unordered_map<int, vector<Point>> &groupedDataP,
                        unordered_map<int, fairofGroup> &fairnessConstraint,  
                        vector<Point> &dataP,                              
                        int groupID,                                        
                        double epsilon,                                    
                        vector<Point> &utiFunClass,                        
                        vector<int> &result,
                        vector<pair<int, double>> &maxUtilityofD, vector<pair<int, double>> &maxUtilityofS,
                        double tauHigh, double tauLow, double &tau, 
                        unordered_map<int, int> &fairofResult,
                        int k, int m, int gamma)
{
      double curDelta = 0.0;
      while (tauHigh - tauLow - epsilon >= 1e-9)
      {
            tau = (tauHigh + tauLow) / 2;
            double sumDelta = 0.0;
            vector<int> tmpResult;
            unordered_map<int, int> tmpFairofResult;
            for (auto it = fairofResult.begin(); it != fairofResult.end(); ++it)
            {
                  tmpFairofResult[it->first] = 0;
            }
            multiroundGreedy(groupedDataP, fairnessConstraint, dataP, groupID, epsilon, utiFunClass, tmpResult,
                             maxUtilityofD, maxUtilityofS, tau, tmpFairofResult, k, m, gamma, sumDelta);
            if (tmpResult.size() != 0 && sumDelta > curDelta)
            {
                  result = tmpResult;
                  fairofResult = tmpFairofResult;
                  curDelta = sumDelta;
                  tauLow = tau;
            }
            else
            {
                  tauHigh = tau;
            }
      }
}

void runBiGreedyAlg(unordered_map<int, vector<Point>> &groupedDataP, 
                      unordered_map<int, fairofGroup> &fairnessConstraint,   
                      vector<Point> &dataP,                               
                      int groupID,                                         
                      double epsilon,                                     
                      vector<Point> &utiFunClass,                          
                      vector<int> &result,
                      int k, int maxM, double &time)
{
      int numofGroups = groupedDataP.size();
      unordered_map<int, int> fairofResult, tmpFairofResult;
      
      for (auto it = fairnessConstraint.begin(); it != fairnessConstraint.end(); ++it)
      {
            fairofResult[it->first] = 0;
      }
      int dim = dataP[0].get_dimension();
      int m; 
      m = k * dim * 5; 
      m = min(m, maxM);
      
      vector<pair<int, double>> maxUtilityofD(m, make_pair(-1, 0.0)), maxUtilityofS(m, make_pair(-1, 0.0)); //int-id double-utility
      for (int i = 0; i < m; ++i)
      {
            double maxUtility = 0.0;
            int maxIdx = -1;
            for (Point p : dataP)
            {
                  int pIdx = p.getId();
                  if (maxUtility < p.dotP(utiFunClass[i]))
                  {
                        maxUtility = p.dotP(utiFunClass[i]);
                        maxIdx = pIdx;
                  }
            }
            maxUtilityofD[i].first = maxIdx;
            maxUtilityofD[i].second = maxUtility;
      }
      
      double tauHigh = 1.0, tauLow = 0.0, tau;

      clock_t start = clock();
      int gamma = log(2 * m / epsilon) / log(2);
      biSearchResult(groupedDataP, fairnessConstraint, dataP, groupID, epsilon, utiFunClass,
                         result, maxUtilityofD, maxUtilityofS, tauHigh, tauLow, tau, fairofResult, k, m, gamma);
                         
      clock_t end = clock();
      time = (double)(end - start) / (CLOCKS_PER_SEC / 1000);
      vector<int> candidate = constructCandidate(groupedDataP, fairnessConstraint, result, fairofResult, k);
      
      while (candidate.size() > 0 && result.size() < k)
      {
            int pIdx = candidate[0];
            Point p = dataP[pIdx];
            result.push_back(pIdx);
            for (int i = 0; i < m; ++i)
            {
                  maxUtilityofS[i].first = p.getId();
                  maxUtilityofS[i].second = max(maxUtilityofS[i].second, p.dotP(utiFunClass[i]));
            }
            ++fairofResult[p.get_category(groupID)];
            candidate = constructCandidate(groupedDataP, fairnessConstraint, result, fairofResult, k);
      }
}

void runBiGreedyPlusAlg(unordered_map<int, vector<Point>> &groupedDataP, 
                          unordered_map<int, fairofGroup> &fairnessConstraint,   
                          vector<Point> &dataP,                               
                          int groupID,                                         
                          double lambda,                                      
                          double epsilon,                                     
                          vector<Point> &utiFunClass,                          
                          vector<int> &result,
                          int k, int maxM, double &time)
{
      int numofGroups = groupedDataP.size();
      unordered_map<int, int> fairofResult; 
      vector<int> tmpResult;
      unordered_map<int, int> tmpFairofResult;
      
	  for (auto it = fairnessConstraint.begin(); it != fairnessConstraint.end(); ++it)
      {
            fairofResult[it->first] = 0;
            tmpFairofResult[it->first] = 0;
      }
      int dim = dataP[0].get_dimension();
      int m; 
      m = dim * k * 1; 
      m = min(m, maxM);
      
      vector<pair<int, double>> maxUtilityofD(m, make_pair(-1, 0.0)), maxUtilityofS(m, make_pair(-1, 0.0)); //int-id double-utility
      for (int i = 0; i < m; ++i)
      {
            double maxUtility = 0.0;
            int maxIdx = -1;
            for (Point p : dataP)
            {
                  int pIdx = p.getId();
                  if (maxUtility < p.dotP(utiFunClass[i]))
                  {
                        maxUtility = p.dotP(utiFunClass[i]);
                        maxIdx = pIdx;
                  }
            }
            maxUtilityofD[i].first = maxIdx;
            maxUtilityofD[i].second = maxUtility;
      }
       
      double tauHigh = 1.0, tauLow = 0.0, tau;
      int gamma = log(2 * m / epsilon) / log(2);
      clock_t start = clock();
      biSearchResult(groupedDataP, fairnessConstraint, dataP, groupID, epsilon, utiFunClass,
                         result, maxUtilityofD, maxUtilityofS, tauHigh, tauLow, tau, fairofResult, k, m, gamma);
      while (true)
      {
            m *= 2;
            m  = min(m, maxM);
            for (int i = m / 2; i < m; ++i)
            {
                  double maxUtility = 0.0;
                  int maxIdx = -1;
                  for (Point p : dataP)
                  {
                        int pIdx = p.getId();
                        if (maxUtility < p.dotP(utiFunClass[i]))
                        {
                              maxUtility = p.dotP(utiFunClass[i]);
                              maxIdx = pIdx;
                        }
                  }
                  maxUtilityofD.push_back(make_pair(maxIdx, maxUtility));
                  maxUtilityofS.push_back(make_pair(-1, 0.0));
            }
            vector<pair<int, double>> tmpMaxS(m, make_pair(-1, 0.0));
            double sumDelta;
            multiroundGreedy(groupedDataP, fairnessConstraint, dataP, groupID, epsilon, utiFunClass, tmpResult,
                             maxUtilityofD, tmpMaxS, tau - lambda, tmpFairofResult, k, m, gamma, sumDelta);
            if (tmpResult.size() != 0)
            {
                  biSearchResult(groupedDataP, fairnessConstraint, dataP, groupID, epsilon, utiFunClass,
                                     result, maxUtilityofD, maxUtilityofS, tau, tau - lambda, tau, fairofResult, k, m, gamma);
                  break;
            }
            else
            {
                  biSearchResult(groupedDataP, fairnessConstraint, dataP, groupID, epsilon, utiFunClass,
                                     result, maxUtilityofD, maxUtilityofS, tau - lambda, 0, tau, fairofResult, k, m, gamma);
            }
      }
      clock_t end = clock();
      time = (double)(end - start) / (CLOCKS_PER_SEC / 1000);
      vector<int> candidate = constructCandidate(groupedDataP, fairnessConstraint, result, fairofResult, k);
      while (candidate.size() > 0 && result.size() < k)
      {
            int pIdx = candidate[0];
            Point p = dataP[pIdx];
            result.push_back(pIdx);
            for (int i = 0; i < m; ++i)
            {
                  maxUtilityofS[i].first = p.getId();
                  maxUtilityofS[i].second = max(maxUtilityofS[i].second, p.dotP(utiFunClass[i]));
            }
            ++fairofResult[p.get_category(groupID)];
            candidate = constructCandidate(groupedDataP, fairnessConstraint, result, fairofResult, k);
      }
}