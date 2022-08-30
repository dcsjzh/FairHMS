#include <bits/stdc++.h>
#include "Point.h"
#include "RMSUtils.h"
#include "greedy.h"
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
    string datasetPath = "../data/", utilsFileName = "../utils/utils_";
    datasetPath += argv[1];
    int k, dim, groupDim, is2D;
    k = atoi(argv[2]);
    is2D = atoi(argv[3]);

    int maxNumofGroups = 5;
    int maxM = 100000;
    vector<unordered_map<string, int>> groupToInt(maxNumofGroups); // transform categorical attribute values to integers, i.e., female -> 0, male -> 1
    vector<Point> dataP;
    readDataPoints(datasetPath, dataP, dim, groupDim, groupToInt);
    vector<Point> utiFunClass;
    utilsFileName += to_string(dim) + "d.txt";
    readUtilityFunctions(utilsFileName, utiFunClass, maxM);


    for (int i = 0; i < groupDim; ++i) 
    {
        double epsilon = 0.02;
        unordered_map<int, vector<Point>> groupedDataP;
        generateGroups(groupedDataP, i, dataP);

        unordered_map<int, fairofGroup> fairnessConstraint = generateFairCons(groupedDataP, dataP, k);
        vector<int> result;
        double time;

        int numofAlgs = 15;
        vector<double> CPUTimes(numofAlgs, -1.0), MHRs(numofAlgs, -1.0);
        double MHR;

        // run the 2D algorithm
        // When the number of groups is greater than 7, IntCov needs too long time to run. We skip it.
        if (is2D == 1 && fairnessConstraint.size() < 7)
        {
            int k = 0;
            vector<int> fairness_upper, fairness_lower;
            double mhr = -1.0;
            for (int groupID = 0; groupID < groupedDataP.size(); groupID++)
            {
                k += fairnessConstraint[groupID].ki;
                fairness_lower.push_back(fairnessConstraint[groupID].lc);
                fairness_upper.push_back(fairnessConstraint[groupID].uc);
            }
            runIntCov(datasetPath, result, {dim + i}, fairness_upper, fairness_lower, k, time, mhr);

            unordered_map<int, int> fairofResult;
            for (auto it = fairnessConstraint.begin(); it != fairnessConstraint.end(); ++it)
            {
                fairofResult[it->first] = 0;
            }
            for (int pIdx : result)
            {
                ++fairofResult[dataP[pIdx].get_category(i)];
            }

            vector<int> candidate = constructCandidate(groupedDataP, fairnessConstraint, result, fairofResult, k);
            while (candidate.size() > 0 && result.size() < k)
            {
                int pIdx = candidate[0];
                Point p = dataP[pIdx];
                result.push_back(pIdx);
                ++fairofResult[p.get_category(i)];
                candidate = constructCandidate(groupedDataP, fairnessConstraint, result, fairofResult, k);
            }
            writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], result, k, i, "IntCov: ", time, mhr);
            MHRs[0] = mhr;
            CPUTimes[0] = time;
        }
        else
        {
            MHRs[0] = -1;
            CPUTimes[0] = -1;
        }
        
        result.clear();
        time = -1.0;
        cout << "Running Bi-Greedy..." << endl;
        runBiGreedyAlg(groupedDataP, fairnessConstraint, dataP, i, epsilon, utiFunClass, result, k, maxM, time);
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], result, k, i, "Bi-Greedy: ", time, MHR);
        MHRs[1] = MHR;
        CPUTimes[1] = time;
        
        result.clear();
        time = -1.0;
        cout << "Running Bi-Greedy-Plus..." << endl;
        runBiGreedyPlusAlg(groupedDataP, fairnessConstraint, dataP, i, 2 * epsilon, epsilon, utiFunClass, result, k, maxM, time);
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], result, k, i, "Bi-Greedy-Plus: ", time, MHR);
        MHRs[2] = MHR;
        CPUTimes[2] = time;

        // run the baseline algorithms

        // run RDP-Greedy on each group separately
        time = 0.0;
        vector<int> resultofGreedy;
        for (int j = 0; j < groupedDataP.size(); ++j)
        {
            vector<Point> curP = groupedDataP[j], curR;
            double timeGreedy = -1.0;
            int r = fairnessConstraint[j].ki;
            runGreedy(curP, r, 1, curR, curP, timeGreedy);
            time += timeGreedy;
            for (Point p : curR)
            {
                resultofGreedy.push_back(p.getId());
            }
        }
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofGreedy, k, i, "Group RDP-Greedy: ", time, MHR);
        MHRs[3] = MHR;
        CPUTimes[3] = time;

        // run greedy directly
        resultofGreedy.clear();
        vector<Point> curR;
        runGreedy(dataP, k, 1, curR, dataP, time);
        for (Point p : curR)
        {
            resultofGreedy.push_back(p.getId());
        }
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofGreedy, k, i, "RDP-Greedy: ", time, MHR);
        MHRs[13] = MHR;
        CPUTimes[13] = time;
        
        // run MRDP-greedy
        resultofGreedy.clear();
        curR.clear();
        runMatroidGreedy(dataP, k, 1, curR, dataP, time, i, groupedDataP, fairnessConstraint);
        for (Point p : curR)
        {
            resultofGreedy.push_back(p.getId());
        }
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofGreedy, k, i, "Fair RDP-Greedy: ", time, MHR);
        MHRs[8] = MHR;
        if(MHR < 0) time = -1;
        CPUTimes[8] = time;


        // run DMM on each group separately
        time = 0.0;
        vector<int> resultofDMMRRMS;
        for (int j = 0; j < groupedDataP.size(); ++j)
        {
            vector<Point> curP = groupedDataP[j], curR;
            double timeDMMRRMS = -1.0;
            int r = fairnessConstraint[j].ki;
            if (r < dim)
            {
                continue;
            }
            if(dim < 10)
                runDMMRRMS(curP, r, 1, curR, curP, timeDMMRRMS);
            time += timeDMMRRMS;
            for (Point p : curR)
            {
                resultofDMMRRMS.push_back(p.getId());
            }
        }
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofDMMRRMS, k, i, "Group DMM-RRMS: ", time, MHR);
        MHRs[4] = MHR;
        CPUTimes[4] = time;
        // run DMM directly
        resultofDMMRRMS.clear();
        curR.clear();
        if(dim < 10)
            runDMMRRMS(dataP, k, 1, curR, dataP, time);
        for (Point p : curR)
        {
            resultofDMMRRMS.push_back(p.getId());
        }
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofDMMRRMS, k, i, "DMM-RRMS: ", time, MHR);
        MHRs[14] = MHR;
        CPUTimes[14] = time;

        // run EPSKernel on each group separately
        time = 0.0;
        vector<int> resultofEps;
        for (int j = 0; j < groupedDataP.size(); ++j)
        {
            vector<Point> curP = groupedDataP[j], curR;
            double timeEps = -1.0;
            int r = fairnessConstraint[j].ki;
            runEpsKernel(curP, r, 1, curR, curP, timeEps);
            time += timeEps;
            for (Point p : curR)
            {
                resultofEps.push_back(p.getId());
            }
        }
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofEps, k, i, "Group Eps-kernel: ", time, MHR);
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
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofEps, k, i, "Eps-kernel: ", time, MHR);
        MHRs[11] = MHR;
        CPUTimes[11] = time;

        // run HS on each group separately
        time = 0.0;
        vector<int> resultofHS;
        for (int j = 0; j < groupedDataP.size(); ++j)
        {
            vector<Point> curP = groupedDataP[j], curR;
            double timeHS = -1.0;
            int r = fairnessConstraint[j].ki;
            runHS(curP, r, 1, curR, curP, timeHS);
            time += timeHS;
            for (Point p : curR)
            {
                resultofHS.push_back(p.getId());
            }
        }
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofHS, k, i, "Group HS: ", time, MHR);
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
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofHS, k, i, "HS: ", time, MHR);
        MHRs[12] = MHR;
        CPUTimes[12] = time;


        // run Sphere on each group separately
        time = 0.0;
        vector<int> resultofSphere;
        for (int j = 0; j < groupedDataP.size(); ++j)
        {
            vector<Point> curP = groupedDataP[j], curR;
            double timeSphere = -1.0;
            int r = fairnessConstraint[j].ki;
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
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofSphere, k, i, "Group Sphere: ", time, MHR);
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
        writeToFile(dataP, dim, groupToInt, fairnessConstraint, argv[1], resultofSphere, k, i, "Sphere: ", time, MHR);
        MHRs[10] = MHR;
        CPUTimes[10] = time;

        MHRs[9] = 0;
        CPUTimes[9] = -1;
        for(int j = 0; j < numofAlgs; j++)
            MHRs[9] = max(MHRs[9], MHRs[j]);

        string outFileName = "../result/";
        outFileName += argv[1];
        outFileName = outFileName.substr(0, outFileName.length() - 4);
        outFileName += "_" + to_string(k) + ".data";
        ofstream fout(outFileName, ios::app);
        if(is2D)
            fout << setw(8) << MHRs[0] << "\t";
        for (int i = 1; i < 10; ++i)
        {
            fout << setw(8) << MHRs[i] << "\t";
        }
        fout << endl;
        if(is2D)
            fout << setw(8) << CPUTimes[0] << "\t";
        for (int i = 1; i < 10; ++i)
        {
            if(MHRs[i] < 0)
                fout << setw(8) << -1 << "\t";
            else
                fout << setw(8) << CPUTimes[i] << "\t";
        }
        fout << endl;
        fout.close();


        if(is2D == 1) continue; 
        // run BiGreedy with ELD

        // Experiments when varying Epsilon and Lambda
        double deltaC = 0.25;
        int setDelta = 1;
        fairnessConstraint = generateFairCons(groupedDataP, dataP, 8);
        vector<vector<double>> CPUTimesEL, MHRsEL;

        int j = 0;
        for(double epsilon = 0.000625; epsilon < 0.7; epsilon *= 2)
        {
            MHRsEL.push_back({});
            CPUTimesEL.push_back({});
            for(double lambda = 0.000625; lambda < 0.7; lambda *= 2)
            {
                
                result.clear();
                time = -1.0;
                cout << "Running Bi-Greedy-Plus..." << endl;
                runBiGreedyPlusAlgWithDelta(groupedDataP, fairnessConstraint, dataP, i, lambda, epsilon, utiFunClass, result, 8, maxM, time, deltaC);
                writeToFileELD(dataP, dim, groupToInt, fairnessConstraint, argv[1], result, 8, i, "Bi-Greedy-Plus: ", time, MHR, epsilon, lambda, deltaC, setDelta);
                MHRsEL[j].push_back(MHR);
                CPUTimesEL[j].push_back(time);
            }
            j++;
        }
        outFileName = "../result/";
        outFileName += argv[1];
        outFileName = outFileName.substr(0, outFileName.length() - 4);
        outFileName += "_" + to_string(8) + "_DC" + to_string(deltaC).substr(0,4) + ".data";
        ofstream foutEL(outFileName, ios::app);
        for(int x = 0; x < MHRsEL.size(); x++)
        {
            for(int y = 0; y < MHRsEL[x].size(); y++)
                foutEL << setw(8) << MHRsEL[x][y] << "\t";
            foutEL << endl;
        }
        for(int x = 0; x < CPUTimesEL.size(); x++)
        {
            for(int y = 0; y < CPUTimesEL[x].size(); y++)
                foutEL << setw(8) << CPUTimesEL[x][y] << "\t";
            foutEL << endl;
        }
        foutEL.close();



        // Experiments when varying Delta
        double lambda = 2*epsilon;
        setDelta = 0;
        fairnessConstraint = generateFairCons(groupedDataP, dataP, 15);
        vector<double> PCPUTimes(numofAlgs, -1.0), PMHRs(numofAlgs, -1.0);

        // run BiGreedy and BiGreedy+ with ELD
        j = 0;
        for(double deltaC = 1.25; deltaC < 60; deltaC *= 2)
        {
            result.clear();
            time = -1.0;
            cout << "Running Bi-Greedy..." << endl;
            runBiGreedyAlgWithDelta(groupedDataP, fairnessConstraint, dataP, i, epsilon, utiFunClass, result, 15, maxM, time, deltaC*10);
            writeToFileELD(dataP, dim, groupToInt, fairnessConstraint, argv[1], result, 15, i, "Bi-Greedy: ", time, MHR, epsilon, lambda, deltaC*10, setDelta);
            MHRs[j] = MHR;
            CPUTimes[j] = time;
            
            result.clear();
            time = -1.0;
            cout << "Running Bi-Greedy-Plus..." << endl;
            runBiGreedyPlusAlgWithDelta(groupedDataP, fairnessConstraint, dataP, i, lambda, epsilon, utiFunClass, result, 15, maxM, time, deltaC);
            writeToFileELD(dataP, dim, groupToInt, fairnessConstraint, argv[1], result, 15, i, "Bi-Greedy-Plus: ", time, MHR, epsilon, lambda, deltaC, setDelta);
            PMHRs[j] = MHR;
            PCPUTimes[j++] = time;
        }

        outFileName = "../result/";
        outFileName += argv[1];
        outFileName = outFileName.substr(0, outFileName.length() - 4);
        outFileName += "_" + to_string(15) + "_E" + to_string(epsilon).substr(0,4) + "_L" + to_string(lambda).substr(0,4) + ".data";
        ofstream foutD(outFileName, ios::app);
        for (int i = 0; i < j; ++i)
            foutD << setw(8) << MHRs[i] << "\t";
        foutD << endl;
        for (int i = 0; i < j; ++i)
            foutD << setw(8) << PMHRs[i] << "\t";
        foutD << endl;
        for (int i = 0; i < j; ++i)
            foutD << setw(8) << CPUTimes[i] << "\t";
        foutD << endl;
        for (int i = 0; i < j; ++i)
            foutD << setw(8) << PCPUTimes[i] << "\t";
        foutD << endl;
        foutD.close();
    }
    return 0;
}