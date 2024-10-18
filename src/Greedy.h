#ifndef GREEDY_H
#define GREEDY_H

#include <iostream>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <optional>
#include "SubsetSelectionAlgorithm.h"
using namespace std;

bool is_streaming = false;

class Greedy : public SubsetSelectionAlgorithm
{
public:
    vector<UtilityFunction> FC;

    Greedy(size_t k, AHR &f, vector<UtilityFunction> &FC) : SubsetSelectionAlgorithm(k, f), FC(FC)
    {

    }

    void fit(vector<Point> const &DataSet, unsigned int iterations = 1)
    {
        vector<size_t> remaining(DataSet.size());
        iota(remaining.begin(), remaining.end(), 0);
        double fcur = 0;                             
        while (solution.size() < k && remaining.size() > 0)
        {
            vector<double> fvals;
            fvals.reserve(remaining.size());
            for (auto i : remaining)
            {
                double ftmp = f.peek(solution, FC, DataSet[i], is_streaming);
                fvals.push_back(ftmp);
            }
            size_t max_fval_pos = distance(fvals.begin(), max_element(fvals.begin(), fvals.end()));
            fcur = fvals[max_fval_pos];
            size_t max_idx = remaining[max_fval_pos];
            f.update(solution, DataSet[max_idx]);
            idxs.push_back(max_idx);
            remaining.erase(remaining.begin() + max_fval_pos);
        }
        fval = fcur;
        is_fitted = true;
    }

    void next(Point const &cur_point, vector<UtilityFunction> &FC)
    {
        throw runtime_error("Please use fit().");
    }


    ~Greedy() {}
};

#endif
