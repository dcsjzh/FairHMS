#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <queue>
#include <map>
#include "Point.h"
#include "RMSUtils.h"
#include "greedyK.h"
#include "DMM.h"
#include "sphere.h"
#include "epskernel.h"
#include "hs.h"
#include "gurobi_c++.h"
#include "IntCov.h"
#include "MyUtils.h"
#include "BiGreedy.h"
using namespace std;

int main(int argc, char *argv[])
{
    string dataset = "../data/", utilsFileName = "../utils/utils_";
    dataset += argv[1];
    int k, dim, categoryDim, is2D;
    categoryDim = atoi(argv[2]);
    k = atoi(argv[3]);
    is2D = atoi(argv[4]);

    int maxNumofCategory = 5;
    int maxM = 10000;
    vector<unordered_map<string, int>> categoryToInt(maxNumofCategory); // transform categorical attribute values to integers, i.e., female -> 0, male -> 1
    vector<Point> dataP;
    readDataPoint(dataset, dataP, dim, categoryDim, categoryToInt);
    vector<Point> utility;
    utilsFileName += to_string(dim) + "d.txt";
    readUtilityFunctions(utilsFileName, utility, maxM);

    double epsilon = 0.02;

    for (int i = 0; i < categoryDim; ++i) 
    {
        unordered_map<int, vector<Point>> categorizedDataP;
        groupGen(categorizedDataP, i, dataP);

        unordered_map<int, fairinCate> partitionMatroid = constructFairCons(categorizedDataP, dataP, k);
        vector<int> result;
        double time;

        int numofAlgs = 13;
        vector<double> CPUTimes(numofAlgs, -1.0), MHRs(numofAlgs, -1.0);
        double MHR;

        // run the 2D algorithm
        if (is2D == 1)
        {
            int k = 0;
            vector<int> fairness_upper, fairness_lower;
            double mhr = -1.0;
            for (int cateId = 0; cateId < categorizedDataP.size(); cateId++)
            {
                k += partitionMatroid[cateId].ki;
                fairness_lower.push_back(partitionMatroid[cateId].lc);
                fairness_upper.push_back(partitionMatroid[cateId].uc);
            }
            if (fairness_lower.size() >= 7)
            {
                continue;
            }
            runIntCov(dataset, result, {dim + i}, fairness_upper, fairness_lower, k, time, mhr);

            unordered_map<int, int> fairofR;
            for (auto it = partitionMatroid.begin(); it != partitionMatroid.end(); ++it)
            {
                fairofR[it->first] = 0;
            }
            for (int pid : result)
            {
                ++fairofR[dataP[pid].get_category(i)];
            }

            vector<int> candidate = constructCandidate(categorizedDataP, partitionMatroid, result, fairofR, k);
            while (candidate.size() > 0 && result.size() < k)
            {
                int id = candidate[0];
                Point p = dataP[id];
                result.push_back(id);
                ++fairofR[p.get_category(i)];
                candidate = constructCandidate(categorizedDataP, partitionMatroid, result, fairofR, k);
            }
            writeToFile2D(dataP, dim, categoryToInt, partitionMatroid, argv[1], result, k, i, "IntCov: ", time, mhr);
            MHRs[0] = mhr;
            CPUTimes[0] = time;

            result.clear();
            time = -1.0;
            cout << "Running Bi-Greedy..." << endl;
            runBiGreedyAlg(categorizedDataP, partitionMatroid, dataP, i, epsilon, utility, result, k, maxM, time);
            writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], result, k, i, "Bi-Greedy: ", time, MHR);
            MHRs[1] = MHR;
            CPUTimes[1] = time;
            
            result.clear();
            time = -1.0;
            cout << "Running Bi-Greedy-Plus..." << endl;
            runBiGreedyPlusAlg(categorizedDataP, partitionMatroid, dataP, i, 2 * epsilon, epsilon, utility, result, k, maxM, time);
            writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], result, k, i, "Bi-Greedy-Plus: ", time, MHR);
            MHRs[2] = MHR;
            CPUTimes[2] = time;
            
        }
        else
        {
            result.clear();
            time = -1.0;
            cout << "Running Bi-Greedy..." << endl;
            runBiGreedyAlg(categorizedDataP, partitionMatroid, dataP, i, epsilon, utility, result, k, maxM, time);
            writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], result, k, i, "Bi-Greedy: ", time, MHR);
            MHRs[1] = MHR;
            CPUTimes[1] = time;
            
            result.clear();
            time = -1.0;
            cout << "Running Bi-Greedy-Plus..." << endl;
            runBiGreedyPlusAlg(categorizedDataP, partitionMatroid, dataP, i, 2 * epsilon, epsilon, utility, result, k, maxM, time);
            writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], result, k, i, "Bi-Greedy-Plus: ", time, MHR);
            MHRs[2] = MHR;
            CPUTimes[2] = time;
        }

        // run the baseline algorithms

        // run RDP-Greedy on each group separately
        time = 0.0;
        vector<int> resultofGreedy;
        for (int j = 0; j < categorizedDataP.size(); ++j)
        {
            vector<Point> curP = categorizedDataP[j], curR;
            double timeGreedy = -1.0;
            int r = partitionMatroid[j].ki;
            runGreedy(curP, r, 1, curR, curP, timeGreedy);
            time += timeGreedy;
            for (Point p : curR)
            {
                resultofGreedy.push_back(p.getId());
            }
        }
        writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], resultofGreedy, k, i, "RDP-Greedy: ", time, MHR);
        MHRs[3] = MHR;
        CPUTimes[3] = time;
        // run greedy
        resultofGreedy.clear();
        vector<Point> curR;
        runGreedy(dataP, k, 1, curR, dataP, time);
        for (Point p : curR)
        {
            resultofGreedy.push_back(p.getId());
        }
        writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], resultofGreedy, k, i, "Total RDP-Greedy: ", time, MHR);
        MHRs[8] = MHR;
        CPUTimes[8] = time;

        // run DMM on each group separately
        time = 0.0;
        vector<int> resultofDMMRRMS;
        for (int j = 0; j < categorizedDataP.size(); ++j)
        {
            vector<Point> curP = categorizedDataP[j], curR;
            double timeDMMRRMS = -1.0;
            int r = partitionMatroid[j].ki;
            if (r < dim)
            {
                continue;
            }
            runDMMRRMS(curP, r, 1, curR, curP, timeDMMRRMS);
            time += timeDMMRRMS;
            for (Point p : curR)
            {
                resultofDMMRRMS.push_back(p.getId());
            }
        }
        writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], resultofDMMRRMS, k, i, "DMM-RRMS: ", time, MHR);
        MHRs[4] = MHR;
        CPUTimes[4] = time;
        // run DMM directly
        resultofDMMRRMS.clear();
        curR.clear();
        runDMMRRMS(dataP, k, 1, curR, dataP, time);
        for (Point p : curR)
        {
            resultofDMMRRMS.push_back(p.getId());
        }
        writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], resultofDMMRRMS, k, i, "Total DMM-RRMS: ", time, MHR);
        MHRs[9] = MHR;
        CPUTimes[9] = time;

        // run EPSKernel on each group separately
        time = 0.0;
        vector<int> resultofEps;
        for (int j = 0; j < categorizedDataP.size(); ++j)
        {
            vector<Point> curP = categorizedDataP[j], curR;
            double timeEps = -1.0;
            int r = partitionMatroid[j].ki;
            runEpsKernel(curP, r, 1, curR, curP, timeEps);
            time += timeEps;
            for (Point p : curR)
            {
                resultofEps.push_back(p.getId());
            }
        }
        writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], resultofEps, k, i, "Eps-kernel: ", time, MHR);
        MHRs[6] = MHR;
        CPUTimes[6] = time;
        // run EPSKernel directly
        resultofEps.clear();
        curR.clear();
        runEpsKernel(dataP, k, 1, curR, dataP, time);
        for (Point p : curR)
        {
            resultofEps.push_back(p.getId());
        }
        writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], resultofEps, k, i, "Total Eps-kernel: ", time, MHR);
        MHRs[11] = MHR;
        CPUTimes[11] = time;

        // run HS on each group separately
        time = 0.0;
        vector<int> resultofHS;
        for (int j = 0; j < categorizedDataP.size(); ++j)
        {
            vector<Point> curP = categorizedDataP[j], curR;
            double timeHS = -1.0;
            int r = partitionMatroid[j].ki;
            runHS(curP, r, 1, curR, curP, timeHS);
            time += timeHS;
            for (Point p : curR)
            {
                resultofHS.push_back(p.getId());
            }
        }
        writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], resultofHS, k, i, "HS: ", time, MHR);
        MHRs[7] = MHR;
        CPUTimes[7] = time;
        // run HS directly
        resultofHS.clear();
        curR.clear();
        runHS(dataP, k, 1, curR, dataP, time);
        for (Point p : curR)
        {
            resultofHS.push_back(p.getId());
        }
        writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], resultofHS, k, i, "Total HS: ", time, MHR);
        MHRs[12] = MHR;
        CPUTimes[12] = time;

        // run Sphere on each group separately
        time = 0.0;
        vector<int> resultofSphere;
        for (int j = 0; j < categorizedDataP.size(); ++j)
        {
            vector<Point> curP = categorizedDataP[j], curR;
            double timeSphere = -1.0;
            int r = partitionMatroid[j].ki;
            if (r < dim)
            {
                continue;
            }
            runSphere(curP, r, 1, curR, curP, timeSphere);
            time += timeSphere;
            for (Point p : curR)
            {
                resultofSphere.push_back(p.getId());
            }
        }
        writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], resultofSphere, k, i, "Sphere: ", time, MHR);
        MHRs[5] = MHR;
        CPUTimes[5] = time;
        // run Sphere directly
        resultofSphere.clear();
        curR.clear();
        runSphere(dataP, k, 1, curR, dataP, time);
        for (Point p : curR)
        {
            resultofSphere.push_back(p.getId());
        }
        writeToFile(dataP, dim, categoryToInt, partitionMatroid, argv[1], resultofSphere, k, i, "Total Sphere: ", time, MHR);
        MHRs[10] = MHR;
        CPUTimes[10] = time;

        string outFileName = "../result/";
        outFileName += argv[1];
        outFileName = outFileName.substr(0, outFileName.length() - 4);
        outFileName += "_" + to_string(k) + ".txt";
        ofstream fout(outFileName, ios::app);
        for (int i = 0; i < numofAlgs; ++i)
        {
            fout << MHRs[i] << "\t";
        }
        fout << endl;
        for (int i = 0; i < numofAlgs; ++i)
        {
            fout << CPUTimes[i] << "\t";
        }
        fout << endl;
        fout.close();
    }
}