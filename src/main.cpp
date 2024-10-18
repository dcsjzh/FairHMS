#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include "Point.h"
#include "IOUtil.h"
#include "RMSUtils.h"
#include "AvgHappiness.h"
#include "Greedy.h"
#include "ThreeSieves.h"
#include "Preemption.h"
#include "GreedyAT.h"
#include "SieveOnline.h"
#include "ThresholdOnline.h"
using namespace std;

const char *FilePath_0 = "./dataset/real/Football_3d/Football_3d_normalized.txt";

ofstream outfile,outfile_fval,outfile_time;

double evaluate_singleton_m(vector<Point> const &DataSet, vector<UtilityFunction> FC, AHR &ahr)
{
    vector<double> ahr_singleton_p(DataSet.size());
    for (size_t i = 0; i < DataSet.size(); ++i)
    {
        double ahr_p = ahr.operator()({DataSet[i]}, FC);
        ahr_singleton_p[i] = ahr_p;
    }
    double max_ahr_singleton_p = *max_element(ahr_singleton_p.begin(), ahr_singleton_p.end());
    return max_ahr_singleton_p;
}

auto evaluate_optimizer(SubsetSelectionAlgorithm &opt, vector<Point> &DataSet)
{
    auto start = chrono::steady_clock::now();
    opt.fit(DataSet);
    auto end = chrono::steady_clock::now();
    chrono::duration<double> runtime_seconds = end - start;
    double fval = opt.get_fval();
    return make_tuple(fval, runtime_seconds.count(),opt.solution.size());
}

int main(int argc, char **argv)
{
    vector<string> FilePaths = {
                                "./dataset/real/Football_3d/Football_3d_normalized.txt",
                                "./dataset/real/Tweet_7d/Tweet_7d_normalized.txt",
                                "./dataset/real/Weather_15d/Weather_15d_normalized.txt",
                                "./dataset/synthetic/IND_3_10000_normalized.txt",
                                // "./dataset/synthetic/Anti-Cor_3_100_normalized.txt",
                                // "./dataset/synthetic/Anti-Cor_3_1000_normalized.txt",
                                "./dataset/synthetic/Anti-Cor_3_10000_normalized.txt",
                                // "./dataset/synthetic/Anti-Cor_3_100000_normalized.txt",
                                // "./dataset/synthetic/Anti-Cor_3_1000000_normalized.txt"
                                // "./dataset/synthetic/Anti-Cor_3_10000_normalized.txt",
                                // "./dataset/synthetic/Anti-Cor_4_10000_normalized.txt",
                                // "./dataset/synthetic/Anti-Cor_5_10000_normalized.txt",
                                // "./dataset/synthetic/Anti-Cor_6_10000_normalized.txt",
                                // "./dataset/synthetic/Anti-Cor_7_10000_normalized.txt",
                                };
    for (size_t count = 0; count < FilePaths.size(); ++count)
    {
        size_t dim;
        vector<Point> dataP;
        IOUtil::read_input_points(FilePaths[count].c_str(), dim, dataP);

        // vector<Point> dataP_temp;
        // dataP_temp.swap(dataP);
        // int offset = 0.2*(count+1)*dataP_temp.size();
        // dataP.assign(dataP_temp.begin(),dataP_temp.begin()+offset);

        /* Read Dataset */
        cout << "Reading data..." << endl;
        cout << "The size of Dataset is " << dataP.size() << endl;
        cout << "The dimension of Dataset is " << dim << endl;

        /* sample ndir utility functions and initialize */
        vector<UtilityFunction> FunctionClass;
        size_t ndir = 1000;
        RMSUtils::get_random_utility_functions(1.0, dim, ndir, FunctionClass, true);
        cout << "The size of FunctionClass is " << FunctionClass.size() << endl;
        for (size_t j = 0; j < FunctionClass.size(); ++j)
        {
            for (size_t i = 0; i < dataP.size(); ++i)
            {
                double f_tmp = FunctionClass[j].direction.dotP(dataP[i]);
                FunctionClass[j].fmax = f_tmp > FunctionClass[j].fmax ? f_tmp : FunctionClass[j].fmax;
            }
        }

        AHR ahr;
        double singleton_m = evaluate_singleton_m(dataP, FunctionClass, ahr); // the largest function value of a singleton set
        tuple<double, double,size_t> res;

        /* the cardinality constraint */
        vector<size_t> ks = {10,20,30,40,50};
        // vector<size_t> ks = {30};

        vector<string> outfiles={
            // "./result/result_for_fig/ahr_varyP_ThresholdOnline_Tweet.txt",
            // "./result/result_for_fig/ahr_varyP_SieveOnline_Tweet.txt"
        };
        vector<string> outfiles_fval={            
            "./result/result_for_fig/Football/ahr_Football_fval.txt",
            "./result/result_for_fig/Tweet/ahr_Tweet_fval.txt",
            "./result/result_for_fig/Weather/ahr_Weather_fval.txt",
            "./result/result_for_fig/IND/ahr_IND_fval.txt",
            "./result/result_for_fig/Anti-Cor/ahr_Anti-Cor_fval.txt",
            // "./result/result_for_fig/ahr_varyD_fval.txt",
            // "./result/result_for_fig/ahr_varyN_fval.txt",
        };
        vector<string> outfiles_time={            
            "./result/result_for_fig/Football/ahr_Football_time.txt",
            "./result/result_for_fig/Tweet/ahr_Tweet_time.txt",
            "./result/result_for_fig/Weather/ahr_Weather_time.txt",
            "./result/result_for_fig/IND/ahr_IND_time.txt",
            "./result/result_for_fig/Anti-Cor/ahr_Anti-Cor_time.txt",
            // "./result/result_for_fig/ahr_varyD_time.txt",
            // "./result/result_for_fig/ahr_varyN_time.txt",
        };
        
        for (size_t count_k = 0; count_k < ks.size(); ++count_k)
        {
            size_t k = ks[count_k];
            // outfile.open(outfiles[count], ios::out|ios::app);
            outfile_fval.open(outfiles_fval[count], ios::out|ios::app);
            outfile_time.open(outfiles_time[count], ios::out|ios::app);
            
            // outfile.open(outfiles[0], ios::out|ios::app);
            // outfile_fval.open(outfiles_fval[0], ios::out|ios::app);
            // outfile_time.open(outfiles_time[0], ios::out|ios::app);
            
            // outfile << "d = " << dim << " n = " << dataP.size() <<  endl;

            // Greedy
            Greedy MyGreedy(k, ahr, FunctionClass);
            res = evaluate_optimizer(MyGreedy, dataP);
            cout << "Selecting " << k << " representative points by Greedy" << "\t fval:\t" << get<0>(res) << "\t runtime:\t" << get<1>(res) << endl;
            // cout << endl;
            // outfile << "Selecting " << k << " representative points by Greedy" << "\t fval:\t" << get<0>(res) << "\t runtime:\t" << get<1>(res) << endl;
            outfile_fval<<get<0>(res)<<" ";
            outfile_time<<get<1>(res)<<" ";            
            // outfile << endl;


            // GreedyAT
            // auto beta_GreedyAT = {0.1, 0.3, 0.5, 0.7, 0.9};
            // auto c1_GreedyAT = {0.6, 0.8, 1.0, 1.2, 1.4};
            // auto c2_GreedyAT = {0.2, 0.4, 0.6, 0.8, 1.0};
            auto beta_GreedyAT = {0.3};
            auto c1_GreedyAT = {1.2};
            auto c2_GreedyAT = {0.6};
            for(auto beta : beta_GreedyAT)
            {
                for(auto c1 : c1_GreedyAT)
                {
                    for(auto c2 : c2_GreedyAT)
                    {
                        if(c2 <= c1)
                        {
                            GreedyAT MyGreedyAT(k, ahr, beta, c1, c2, FunctionClass);
                            res = evaluate_optimizer(MyGreedyAT, dataP);
                            cout << "Selecting " << k << " representatives by GreedyAT with beta = " << beta << " and c1 = " << c1 << " and c2 = " << c2 << "\t fval:\t" << get<0>(res) << "\t runtime:\t" << get<1>(res) << endl;
                            outfile_fval<<get<0>(res)<<" ";
                            outfile_time<<get<1>(res)<<" "; 
                            // outfile << "Selecting " << k << " representatives by GreedyAT with beta = " << beta << " and c1 = " << c1 << " and c2 = " << c2 << "\t fval:\t" << get<0>(res) << "\t runtime:\t" << get<1>(res) << endl;
                        }
                    }
                }
            }


            // Preemption
            // auto c_Preemption = {0.6, 0.8, 1.0, 1.2, 1.4};
            auto c_Preemption = {1.0};
            for(auto c: c_Preemption)
            {
                Preemption MyPreemption(k, ahr, c, FunctionClass);
                res = evaluate_optimizer(MyPreemption, dataP);
                cout << "Selecting " << k << " representatives by Preemption with c = " << c << "\t fval:\t" << get<0>(res) << "\t runtime:\t" << get<1>(res) << endl;
                // outfile << "Preemption Selecting " << k << "\t fval:\t" << get<0>(res) << "\t runtime:\t" << get<1>(res) << endl;
                outfile_fval<<get<0>(res)<<" ";
                outfile_time<<get<1>(res)<<" ";    
            }


            // ThreeSieves
            // auto T_ThreeSieves = {50, 100, 500, 1000, 5000};
            // auto eps_ThreeSieves = {0.001, 0.005, 0.01, 0.05, 0.1};
            auto T_ThreeSieves = {1000};
            auto eps_ThreeSieves = {0.1};
            for(auto T : T_ThreeSieves)
            {
                for(auto eps : eps_ThreeSieves)
                {
                    ThreeSieves MyThreeSieves(k, ahr, singleton_m, eps, ThreeSieves::THRESHOLD_STRATEGY::SIEVE, T, FunctionClass);
                    res = evaluate_optimizer(MyThreeSieves, dataP);
                    cout << "Selecting " << k << " representatives by ThreeSieves with T = " << T << " and epsilon = " << eps  << "\t fval:\t" << get<0>(res) << "\t runtime:\t" << get<1>(res) << endl;
                    // outfile << "Selecting " << k << " representatives by ThreeSieves with T = " << T << " and epsilon = " << eps  << "\t fval:\t" << get<0>(res) << "\t runtime:\t" << get<1>(res) << endl;
                    outfile_fval<<get<0>(res)<<" ";
                    outfile_time<<get<1>(res)<<" "; 
                }
            }


            // ThresholdOnline
            // auto alpha_ThresholdOnline = {0.0, 0.02, 0.04, 0.06, 0.08, 0.1};
            // auto gamma_ThresholdOnline = {0.0,0.25,0.5, 0.75,1.0, 1.25, 1.5};
            auto alpha_ThresholdOnline = {0.0};
            auto gamma_ThresholdOnline = {0.0};
            for (auto alpha : alpha_ThresholdOnline)
            {
                for(auto gamma : gamma_ThresholdOnline)
                {
                    ThresholdOnline MyThresholdOnline(k, ahr, alpha, gamma, FunctionClass);
                    res = evaluate_optimizer(MyThresholdOnline, dataP);
                    cout << "Selecting " << k <<"->"<<get<2>(res)<< " representatives by ThresholdOnline with alpha = " << alpha << " and gamma = " << gamma << "\t fval:\t" << get<0>(res) << "\t runtime:\t" << get<1>(res) << endl;
                    // outfile <<  alpha << " " << gamma << " "<< get<0>(res) << " "<< get<1>(res) << " "<< get<2>(res)<< endl;
                    outfile_fval<<get<0>(res)<<" ";
                    outfile_time<<get<1>(res)<<" "; 
                }
            }


            // SieveOnline
            // auto eps_SieveOnline = {0.002, 0.005, 0.02, 0.05, 0.2};
            auto eps_SieveOnline = {0.2};
            for (auto eps : eps_SieveOnline)
            {
                SieveOnline MySieveOnline(k, ahr, eps, FunctionClass);
                res = evaluate_optimizer(MySieveOnline, dataP);
                cout << "Selecting " << k <<"->"<<get<2>(res) << " representatives by SieveOnline with eps = " << eps << "\t fval:\t" << get<0>(res) << "\t\t runtime:\t" << get<1>(res) << endl;
                // outfile << "SieveOnline Selecting " << k <<"->"<<get<2>(res) << "\t fval:\t" << get<0>(res) << "\t\t runtime:\t" << get<1>(res) << endl;
                // outfile << eps << " " << get<0>(res) << " " << get<1>(res) << endl;
                outfile_fval<<get<0>(res)<<" ";
                outfile_time<<get<1>(res)<<" "; 
            }

            cout << endl;
            // outfile << endl;
            outfile_fval<<endl;
            outfile_time<<endl; 
            // outfile.close();
            outfile_fval.close();
            outfile_time.close();
        }
    }
    return 0;
}