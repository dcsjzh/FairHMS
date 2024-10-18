#ifndef THRESHOLDONLINE_H
#define THRESHOLDONLINE_H

#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
#include <unordered_set>
#include "SubsetSelectionAlgorithm.h"
using namespace std;

class ThresholdOnline : public SubsetSelectionAlgorithm
{
public:
    double alpha;
    double beta;
    vector<UtilityFunction> FC;
    bool is_streaming = true;

    ThresholdOnline(size_t k, AHR &f, double alpha, double beta, vector<UtilityFunction> &FC) : SubsetSelectionAlgorithm(k, f), alpha(alpha), beta(beta), FC(FC)
    {
        
    }

    void fit(vector<Point> const &DataSet, unsigned int iterations = 1)
    {
        size_t N = DataSet.size();
        for (size_t t = 0; t < N; ++t)
        {
            double fdelta = f.peek(solution, FC, DataSet[t], is_streaming) - fval;
            if(fdelta > alpha + (beta + 1) * fval / k)
            {
                size_t cur_size = solution.size();
                if(cur_size == k)
                {
                    vector<double> fvals;       
                    fvals.reserve(cur_size);
                    vector<size_t> idxs_remaining(idxs);
                    vector<Point> set_tmp;
                    double f_pre = 0.0, f_last;
                    for(size_t j = 0; j < k; ++j)
                    {
                        size_t min_idxs_pos = distance(idxs_remaining.begin(), min_element(idxs_remaining.begin(), idxs_remaining.end()));
                        f_last = f.peek(set_tmp, FC, solution[min_idxs_pos], is_streaming);
                        fvals[min_idxs_pos] = f_last - f_pre;
                        set_tmp.push_back(solution[min_idxs_pos]);
                        idxs_remaining.erase(idxs_remaining.begin() + min_idxs_pos);
                        f_pre = f_last;
                    }
                    size_t min_fval_pos = distance(fvals.begin(), min_element(fvals.begin(), fvals.end()));
                    solution[min_fval_pos] = DataSet[t];
                    idxs[min_fval_pos] = t;
                    fval = f.operator()(solution, FC);
                }
                else
                {
                    f.update(solution, DataSet[t]);
                    fval += fdelta;
                    idxs.push_back(t);
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