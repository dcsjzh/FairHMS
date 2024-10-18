#ifndef SIEVEONLINE_H
#define SIEVEONLINE_H

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <numeric>
#include <unordered_set>
using namespace std;

inline vector<double> thresholds(double lower, double upper, double epsilon)
{
    vector<double> ts;
    if (epsilon > 0.0)
    {
        int ilower = ceil(log(lower) / log(1.0 + epsilon));
        for (double val = pow(1.0 + epsilon, ilower); val <= upper; ++ilower, val = pow(1.0 + epsilon, ilower))
        {
            ts.push_back(val);
        }
    }
    else
    {
        throw runtime_error("thresholds: epsilon must be a positive real-number (is: " + to_string(epsilon) + ").");
    }
    return ts;
}

class SieveOnline : public SubsetSelectionAlgorithm
{
private:
    class Sieve : public SubsetSelectionAlgorithm
    {
    public:
        double threshold;
        bool is_streaming = true;
        size_t Counter = 0;

        Sieve(size_t k, AHR &f, double threshold) : SubsetSelectionAlgorithm(k, f), threshold(threshold) 
        {

        }

        void fit(vector<Point> const &DataSet, unsigned int iterations = 1)
        {
            throw runtime_error("Sieves are only meant to be used through SieveOnline and therefore do not require the implementation of `fit'");
        }

        void next(Point const &cur_point, vector<UtilityFunction> &FC)
        {
            size_t cur_size = solution.size();
            if (cur_size < k)
            {
                double fdelta = f.peek(solution, FC, cur_point, is_streaming) - fval;
                double tau = (threshold / 2.0 - fval) / static_cast<double>(k - cur_size);
                if (fdelta >= tau)
                {
                    f.update(solution, cur_point);
                    fval += fdelta;
                }
            }
            is_fitted = true;
        }
    };

public:
    list<unique_ptr<Sieve>> sieves;
    vector<UtilityFunction> FC;
    bool is_streaming = true;
    double epsilon;
    double m = 0;
    SieveOnline(size_t k, AHR &f, double epsilon, vector<UtilityFunction> &FC) : SubsetSelectionAlgorithm(k, f), FC(FC), epsilon(epsilon)
    {

    }

    void fit(vector<Point> const &DataSet, unsigned int iterations = 1)
    {
        for (auto &cur_point : DataSet)
        {
            next(cur_point, FC);
        }
    }

    void next(Point const &cur_point, vector<UtilityFunction> &FC)
    {
        bool m_isChange = false;
        double m_temp = f.operator()({cur_point}, FC);
        if(m_temp > m)
        {
            m = m_temp;
            m_isChange = true;
        }
        if(m_isChange)
        {
            vector<double> ts = thresholds(m, 2*k*m, epsilon);
            while(!sieves.empty() && sieves.front()->threshold != ts.front())
            {
                sieves.erase(sieves.begin());
            }
            list<unique_ptr<Sieve>>::iterator it_sieves = sieves.begin();
            vector<double>::iterator it_ts = ts.begin();
            while(it_sieves != sieves.end())
            {
                ++it_sieves;
                ++it_ts;
            }
            while(it_ts != ts.end())
            {
                sieves.push_back(make_unique<Sieve>(k, f, *it_ts));
                ++it_ts;
            }          
        }
        bool fval_sieves_isChange = false;
        for (auto &s : sieves)
        {
            s->next(cur_point, FC);
            if (s->get_fval() > fval_sieves)
            {
                fval_sieves = s->get_fval();
                fval_sieves_isChange = true;
            }
        }  
        size_t cur_size = solution.size();
        if(cur_size < k)
        {
            double fdelta = f.peek(solution, FC, cur_point, is_streaming) - fval;
            if (fdelta >= fval_sieves / (2 * k))
            {
                f.update(solution, cur_point);
                fval += fdelta;
            }
        }
        else
        {
            vector<Point> tmp_solution(solution);
            tmp_solution.push_back(cur_point);
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
                if (l == tmp_solution.size()-1)
                {
                    fval = sum / FC.size();  
                }
                else   
                {
                    double ftmp = sum / FC.size();
                    fvals.push_back(ftmp);
                }
            }
            size_t max_fval_pos = distance(fvals.begin(), max_element(fvals.begin(), fvals.end()));
            double fcur = fvals[max_fval_pos];
            if (fcur - fval >= fval_sieves / (2 * k))
            {
                solution[max_fval_pos] = cur_point;
                fval = fcur;
            }
        }
        is_fitted = true;
    };

    unsigned int get_num_candidate_solutions() const
    {
        return sieves.size();
    }

    unsigned long get_num_elements_stored() const
    {
        unsigned long num_elements = 0;
        for (auto const &s : sieves)
        {
            num_elements += s->get_solution().size();
        }

        return num_elements;
    }

    ~SieveOnline() {}
};

#endif