#ifndef THREESIEVES_H
#define THREESIEVES_H

#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
#include <unordered_set>
#include <string>
#include "SubsetSelectionAlgorithm.h"
using namespace std;

class ThreeSieves : public SubsetSelectionAlgorithm
{
public:
    enum THRESHOLD_STRATEGY
    {
        SIEVE,   /*!< Start with the largest threshold in \f$ \{(1+\varepsilon)^i  | i \in Z, lower \le (1+\varepsilon)^i \le upper\} \f$ and always use the next largest as the new threshold */
        CONSTANT /*!< Reduce the threshold by a constant \f$ \varepsilon \f$  */
    };

    double threshold;
    double epsilon;
    THRESHOLD_STRATEGY strategy;
    unsigned int T;
    unsigned int t;
    vector<UtilityFunction> FC;
    bool is_streaming = true;

    ThreeSieves(size_t k, AHR &f, double m, double epsilon, THRESHOLD_STRATEGY strategy, size_t T, vector<UtilityFunction> &FC) : SubsetSelectionAlgorithm(k, f), threshold(k * m), epsilon(epsilon), strategy(strategy), T(T), t(0), FC(FC)
    {
        assert(T >= 1);
    }

    void fit(vector<Point> const &DataSet, unsigned int iterations = 1)
    {
        for (auto &cur_point : DataSet)
        {
            next(cur_point, FC);
            if (solution.size() == k)
            {
                return;
            }
        }
    }

    void next(Point const &cur_point, vector<UtilityFunction> &FC)
    {
        size_t cur_size = solution.size();
        if (cur_size < k)
        {
            if (t >= T)
            {
                switch (strategy)
                {
                case THRESHOLD_STRATEGY::SIEVE:
                {
                    double tmp = log(threshold) / log(1.0 + epsilon);
                    int i;
                    if (tmp == floor(tmp) || abs(tmp - floor(tmp)) < 1e-7)
                    {
                        i = floor(tmp) - 1;
                    }
                    else
                    {
                        i = floor(tmp);
                    }
                    threshold = pow(1 + epsilon, i);
                    break;
                }
                case THRESHOLD_STRATEGY::CONSTANT:
                {
                    threshold = threshold - epsilon;
                    break;
                }
                }
                t = 0;
            }
            double fdelta = f.peek(solution, FC, cur_point, is_streaming) - fval;
            double tau = (threshold / 2.0 - fval) / static_cast<double>(k - cur_size);
            if (fdelta >= tau)
            {
                f.update(solution, cur_point);
                fval += fdelta;
                t = 0;
            }
            else
            {
                ++t;
            }
        }
        is_fitted = true;
    }
};

#endif