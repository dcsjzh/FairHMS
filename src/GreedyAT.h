#ifndef GREEDYAT_H
#define GREEDYAT_H

#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
#include <unordered_set>
#include "SubsetSelectionAlgorithm.h"
using namespace std;

class GreedyAT : public SubsetSelectionAlgorithm
{
public:
    double beta;
    double c1;
    double c2;
    vector<UtilityFunction> FC;
    bool is_streaming = true;

    GreedyAT(size_t k, AHR &f, double beta, double c1, double c2, vector<UtilityFunction> &FC) : SubsetSelectionAlgorithm(k, f), beta(beta), c1(c1), c2(c2), FC(FC)
    {
        
    }

    void fit(vector<Point> const &DataSet, unsigned int iterations = 1)
    {
        size_t N = DataSet.size();
        for (size_t t = 0; t < N; ++t)
        {
            if (t <= k - 1)
            {
                f.update(solution, DataSet[t]); 
                if(t == k - 1)
                {
                    fval = f.operator()(solution, FC);
                }
            }
            else
            {
                double c;
                if (t <= beta * N)
                {
                    c = c1;
                }
                else
                {
                    c = c2;
                }
                vector<Point> tmp_solution(solution);
                tmp_solution.push_back(DataSet[t]); 
                vector<vector<double>> happinesses;
                for (size_t j = 0; j < FC.size(); ++j)
                {
                    vector<double> happiness;
                    double happy_temp;
                    for (size_t i = 0; i < tmp_solution.size(); ++i)
                    {
                        happy_temp = FC[j].direction.dotP(tmp_solution[i]);
                        happiness.push_back(happy_temp); 
                    }
                    happinesses.push_back(happiness);
                }
                vector<double> fvals;
                fvals.reserve(solution.size());
                for (size_t l = 0; l < tmp_solution.size(); ++l)
                {
                    double happy_max, happy_ratio;
                    double sum = 0;
                    for (size_t j = 0; j < happinesses.size(); ++j)
                    {
                        happy_max = find_max_happy(happinesses[j], l);
                        happy_ratio = happy_max / FC[j].fmax;
                        sum += happy_ratio;
                    }
                    if (l < tmp_solution.size()-1)
                    {
                        double ftmp = sum / FC.size();
                        fvals.push_back(ftmp);
                    }
                }
                size_t max_fval_pos = distance(fvals.begin(), max_element(fvals.begin(), fvals.end()));
                double fcur = fvals[max_fval_pos];
                if (fcur - fval >= c * fval / k)
                {
                    solution[max_fval_pos] = DataSet[t]; 
                    fval = fcur;
                }
            }
        }
        is_fitted = true;
    }

    void next(Point const &cur_point, vector<UtilityFunction> &FC)
    {
        throw runtime_error("Please use fit().");
    }
};

#endif