#ifndef AVERAGE_HAPPINESS_RATIO_H
#define AVERAGE_HAPPINESS_RATIO_H

#include <mutex>
#include <vector>
#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>
#include "ObjectiveFunction.h"
#include "Point.h"
#include "RMSUtils.h"
using namespace std;

class AHR
{
public:
    double ahr_val;

    AHR() : ahr_val(0) {}

    double operator()(vector<Point> const &cur_solution, vector<UtilityFunction> const &FunctionClass)
    {
        double sum = 0;
        double ratio = 0;
        for (size_t j = 0; j < FunctionClass.size(); ++j)
        {
            double happy_max = 0;
            double happy_temp = 0;
            for (size_t i = 0; i < cur_solution.size(); ++i)
            {
                happy_temp = cur_solution[i].dotP(FunctionClass[j].direction);
                if (happy_temp > happy_max)
                {
                    happy_max = happy_temp;
                }
            }
            ratio = happy_max / FunctionClass[j].fmax;
            sum += ratio;
        }
        ahr_val = sum / FunctionClass.size();
        return ahr_val;
    }

    double peek(vector<Point> const &cur_solution, vector<UtilityFunction> const &FunctionClass, Point const &cur_point, bool is_streaming)
    {
        vector<Point> tmp_solution(cur_solution);
        tmp_solution.push_back(cur_point);
        double ftemp;
        ftemp = this->operator()(tmp_solution, FunctionClass);
        return ftemp;
    }

    void update(vector<Point> &cur_solution, Point const &cur_point)
    {
        cur_solution.push_back(cur_point);
    }

    double get_ahr_val()
    {
        return ahr_val;
    }

    ~AHR() {}
};

#endif