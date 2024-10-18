#ifndef SUBMODULAROPTIMIZER_H
#define SUBMODULAROPTIMIZER_H

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include <functional>
#include <cassert>
#include <optional>
#include <iostream>
#include "ObjectiveFunction.h"
#include "AvgHappiness.h"
using namespace std;

class SubsetSelectionAlgorithm
{
private:
protected:
    size_t k;
    AHR f;
    bool is_fitted;

public:
    vector<Point> solution;
    double fval;
    double fval_sieves;
    vector<size_t> idxs;

    SubsetSelectionAlgorithm(size_t k, AHR &f) : k(k), f(f)
    {
        is_fitted = false;
        fval = 0;
        fval_sieves = 0;
    }

    virtual void fit(vector<Point> const &DataSet, unsigned int iterations = 1)
    {
        for (unsigned int i = 1; i <= iterations; ++i)
        {
            for (auto &cur_point : DataSet)
            {
                vector<UtilityFunction> FC_virtual;
                next(cur_point, FC_virtual);
                if (solution.size() == k && i > 1)
                {
                    return;
                }
            }
        }
    }

    virtual void next(Point const &cur_point, vector<UtilityFunction> &FC) = 0;

    vector<Point> const &get_solution() const
    {
        if (!this->is_fitted)
        {
            throw runtime_error("Optimizer was not fitted yet! Please call fit() or next() before calling get_solution()");
        }
        else
        {
            return solution;
        }
    }

    virtual unsigned int get_num_candidate_solutions() const
    {
        return 1;
    }

    virtual unsigned long get_num_elements_stored() const
    {
        return this->get_solution().size();
    }

    double get_fval() const
    {
        return fval;
    }

    void Update_FC(vector<Point> const &DataSet, Point const &cur_point, vector<UtilityFunction> &FC, bool is_streaming)
    {
        if (!is_streaming)
        {
            for (size_t j = 0; j < FC.size(); ++j)
            {
                FC[j].fmax = 0;
                for (size_t i = 0; i < DataSet.size(); ++i)
                {
                    FC[j].fmax = FC[j].direction.dotP(DataSet[i]) > FC[j].fmax ? FC[j].direction.dotP(DataSet[i]) : FC[j].fmax;
                }
            }
        }
    }

    double find_max_happy(vector<double> &happiness, size_t l)
    {
        double happy_max = 0;

        for (size_t i = 0; i < happiness.size(); ++i)
        {
            if (i == l)
            {
                continue;
            }
            else
            {
                happy_max = happiness[i] > happy_max ? happiness[i] : happy_max;
            }
        }
        return happy_max;
    }

    vector<Point> extract(vector<Point> solution, size_t it)
    {
        vector<Point> tmp_set;
        for(size_t i = 0; i < it; ++i)
        {
            tmp_set.push_back(solution[i]);
        }
        return tmp_set;
    }

    virtual ~SubsetSelectionAlgorithm() {}
};

#endif
