#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include <functional>
#include <cassert>
#include <iostream>
using namespace std;

class ObjectiveFunction {
public:
    virtual double operator()(vector<vector<double>> const &cur_solution) const = 0;

    virtual double peek(vector<vector<double>> const &cur_solution, vector<double> const &x, unsigned int pos) = 0; 

    virtual void update(vector<vector<double>> const &cur_solution, vector<double> const &x, unsigned int pos) = 0;

    virtual shared_ptr<ObjectiveFunction> clone() const = 0;

    virtual ~ObjectiveFunction() {}
};

class ObjectiveFunctionWrapper : public ObjectiveFunction {
protected:
    function<double (vector<vector<double>> const &)> f;

public:
    ObjectiveFunctionWrapper(function<double (vector<vector<double>> const &)> f) : f(f) {}

    double operator()(vector<vector<double>> const &cur_solution) const 
    {
        return f(cur_solution);
    }

    double peek(vector<vector<double>> const &cur_solution, vector<double> const &x, unsigned int pos)
    {
        vector<vector<double>> tmp(cur_solution);
        if (pos >= cur_solution.size())
        {
            tmp.push_back(x);
        }
        else
        {
            tmp[pos] = x;
        }
        double ftmp = this->operator()(tmp);
        return ftmp;
    }

    void update(vector<vector<double>> const &cur_solution, vector<double> const &x, unsigned int pos) {}
    
    shared_ptr<ObjectiveFunction> clone() const 
    {
        return shared_ptr<ObjectiveFunction>(new ObjectiveFunctionWrapper(f));
    }

    ~ObjectiveFunctionWrapper() {}
};

#endif
