#include "MyUtils.h"

vector<string> split(string src, char ch)
{
    int i = 0;
    vector<string> splitString;
    string tmpString;
    while (i < src.length())
    {
        if (src[i] != ch)
        {
            tmpString += src[i++];
        }
        else
        {
            splitString.push_back(tmpString);
            tmpString = "";
            ++i;
        }
    }
    splitString.push_back(tmpString);
    return splitString;
}

void readDataPoints(string fileName, vector<Point> &dataP, int &dim, int &groupDim,
                   vector<unordered_map<string, int>> &groupToInt)
{
    ifstream fin(fileName);
    if (fin)
    {
        string line;
        // the first line represents the dimensionality
        getline(fin, line);
        dim = stoi(line);
        int maxNumofGroups = 5;
        vector<int> groupIDs(maxNumofGroups, 0); 
        int count = 0;
        while (getline(fin, line))
        {
            vector<double> coords;
            vector<int> groups;
            vector<string> tmpString = split(line, '\t');
            groupDim = tmpString.size() - dim;
            for (int j = 0; j < dim; ++j)
            {
                coords.push_back(stod(tmpString[j]));
            }
            for (int j = 0; j < groupDim; ++j)
            {
                if (groupToInt[j].count(tmpString[dim + j]) == 0)
                {
                    groupToInt[j][tmpString[dim + j]] = groupIDs[j]++;
                }
                groups.push_back(groupToInt[j][tmpString[dim + j]]);
            }
            Point p(dim, count++, coords, groups);
            dataP.push_back(p);
        }
    }
    else
    {
        cout << "open \"" + fileName + "\" failed." << endl;
    }
}

/* read the utility functions from the file */
void readUtilityFunctions(string fileName, vector<Point> &utiFunClass, int m)
{
    ifstream fin(fileName);
    if (fin)
    {
        string line;
        int count = 0;
        while (getline(fin, line))
        {
            vector<double> coords;
            for (string s : split(line, ' '))
            {
                coords.push_back(stod(s));
            }
            Point utility(coords.size(), count++, coords);
            utiFunClass.push_back(utility);
            if (count == m)
            {
                break;
            }
        }
    }
}

/* generate groups */
void generateGroups(unordered_map<int, vector<Point>> &groupedDataP,
                   int groupID, vector<Point> dataP)
{
    for (Point p : dataP)
    {
        int pGroupID = p.get_category(groupID);
        if (groupedDataP.count(pGroupID) == 0)
        {
            vector<Point> tmp;
            tmp.push_back(p);
            groupedDataP[pGroupID] = tmp;
        }
        else
        {
            groupedDataP[pGroupID].push_back(p);
        }
    }
}

/* write the results to the folder “result” */
void writeToFile(vector<Point> &dataP, int dim, vector<unordered_map<string, int>> &groupToInt, unordered_map<int, fairofGroup> &fairnessConstraint,
                 string datasetPath, vector<int> &result, int k, int groupID, string algName, double time, double &MHR)
{
    // calculate the MHR of the final solution
    size_t ndir = RMSUtils::ndir_for_validation(dim);
    vector<Point> resultP;
    for (int pIdx : result)
    {
        resultP.push_back(dataP[pIdx]);
    }
    double MaxR = 1.0, AvgR, perc80;
    cout << "Computing HR ..." << endl;
    if (resultP.size() == 0)
    {
        MaxR = 0.0;
        MHR = - 1.0;
    }
    else
    {
        RMSUtils::Max_Avg_Regret(1.0, dim, ndir, dataP, resultP, MaxR, AvgR, perc80, 1);

        MHR = 1.0 - MaxR;
    }

    string outFileName = "../result/";
    outFileName += datasetPath;
    outFileName = outFileName.substr(0, outFileName.length() - 4);
    outFileName += "_" + to_string(k) + ".txt";
    ofstream fout(outFileName, ios::app);
    fout << algName << "groupID=" << groupID << endl
         << "Fairness constraint: ";
    
    vector<string> groupNames;
    for (int j = 0; j < fairnessConstraint.size(); ++j)
    {
        for (auto it = groupToInt[groupID].begin(); it != groupToInt[groupID].end(); ++it)
        {
            if (it->second == j)
            {
                groupNames.push_back(it->first);
                break;
            }
        }
    }
    
    for (int i = 0; i < fairnessConstraint.size(); ++i)
    {
        fout << groupNames[i] << "=" << fairnessConstraint[i].lc << "-" << fairnessConstraint[i].ki << "-" << fairnessConstraint[i].uc << " ";
    }
    fout << endl;
    vector<int> numofCates(fairnessConstraint.size(), 0);
    for (int j = 0; j < fairnessConstraint.size(); ++j)
    {
        string groupNameStr = groupNames[j];
        for (int pIdx : result)
        {
            if (dataP[pIdx].get_category(groupID) == j)
            {
                fout << "idx=" << pIdx << "\t" << groupNameStr << "\t";
                for (int k = 0; k < dim; ++k)
                {
                    fout << to_string(dataP[pIdx].get_coordinate(k)) << "\t";
                }
                fout << endl;
                ++numofCates[j];
            }
        }
    }
    
    for (int i = 0; i < fairnessConstraint.size(); ++i)
    {
        fout << groupNames[i] << "=" << numofCates[i] << " ";
    }
    fout << endl;
    fout << "HR=" << 1.0 - MaxR << endl
         << "Time=" << time << endl
         << endl;
    fout.close();
}

/* write the results of the 2D algorithm - IntCov */
void writeToFile2D(vector<Point> &dataP, int dim, vector<unordered_map<string, int>> &groupToInt, unordered_map<int, fairofGroup> &fairnessConstraint,
                   string datasetPath, vector<int> &result, int k, int groupID, string algName, double time, double mhr)
{
    string outFileName = "../result/";
    outFileName += datasetPath;
    outFileName = outFileName.substr(0, outFileName.length() - 4);
    outFileName += "_" + to_string(k) + ".txt";
    ofstream fout(outFileName, ios::app);
    fout << algName << "groupID=" << groupID << endl
         << "Fairness constraint: ";
    
    vector<string> groupNames;
    for (int j = 0; j < fairnessConstraint.size(); ++j)
    {
        for (auto it = groupToInt[groupID].begin(); it != groupToInt[groupID].end(); ++it)
        {
            if (it->second == j)
            {
                groupNames.push_back(it->first);
                break;
            }
        }
    }

    for (int i = 0; i < fairnessConstraint.size(); ++i)
    {
        fout << groupNames[i] << "=" << fairnessConstraint[i].lc << "-" << fairnessConstraint[i].ki << "-" << fairnessConstraint[i].uc << " ";
    }
    fout << endl;
    vector<int> numofCates(fairnessConstraint.size(), 0);
    for (int j = 0; j < fairnessConstraint.size(); ++j)
    {
        string groupNameStr = groupNames[j];
        for (int pIdx : result)
        {
            if (dataP[pIdx].get_category(groupID) == j)
            {
                fout << "idx=" << pIdx << "\t" << groupNameStr << "\t";
                for (int k = 0; k < dim; ++k)
                {
                    fout << to_string(dataP[pIdx].get_coordinate(k)) << "\t";
                }
                fout << endl;
                ++numofCates[j];
            }
        }
    }
 
    for (int i = 0; i < fairnessConstraint.size(); ++i)
    {
        fout << groupNames[i] << "=" << numofCates[i] << " ";
    }
    fout << endl;
    fout << "HR=" << mhr << endl
         << "Time=" << time << endl
         << endl;
    fout.close();
}

/* generate the fairness constraints */
unordered_map<int, fairofGroup> generateFairCons(unordered_map<int, vector<Point>> groupedDataP, vector<Point> dataP, int k)
{
    unordered_map<int, fairofGroup> fairnessConstraint;
    
    vector<double> propotion(groupedDataP.size());
    int total = 0, smallestGroupID, smallestGroupSize = k + 1;
    for (int i = 0; i < propotion.size(); ++i)
    {
        propotion[i] = groupedDataP[i].size() / (dataP.size() * 1.0);
        fairnessConstraint[i].ki = round(k * propotion[i]);
        total += fairnessConstraint[i].ki;
        if (fairnessConstraint[i].ki < smallestGroupSize)
        {
            smallestGroupSize = fairnessConstraint[i].ki;
            smallestGroupID = i;
        }
    }
    while (total < k)
    {
        ++fairnessConstraint[smallestGroupID].ki;
        ++total;
        smallestGroupSize = fairnessConstraint[smallestGroupID].ki;
        //
        for (int i = 0; i < propotion.size(); ++i)
        {
            if (fairnessConstraint[i].ki < smallestGroupSize)
            {
                smallestGroupSize = fairnessConstraint[i].ki;
                smallestGroupID = i;
            }
        }
    }
    
    int count = 0;
    for (auto it = fairnessConstraint.begin(); it != fairnessConstraint.end(); ++it)
    {
        if (it->second.ki == 0)
        {
            ++(it->second.ki);
            ++count;
        }
    }
    while (count != 0)
    {
        int largestGroupID, largestGroupSize = 0;
        for (auto it = fairnessConstraint.begin(); it != fairnessConstraint.end(); ++it)
        {
            if (it->second.ki > largestGroupSize)
            {
                largestGroupID = it->first;
                largestGroupSize = it->second.ki;
            }
        }
        if (largestGroupSize > 1)
        {
            --fairnessConstraint[largestGroupID].ki;
            --count;
        }
        else
        {
            break;
        }
    }

    int n = dataP.size();
    double alpha = 0.1;
    int gap = alpha * k;
    gap = max(1, gap);
    // set l_c and h_c
    for (auto it = fairnessConstraint.begin(); it != fairnessConstraint.end(); ++it)
    {
        it->second.lc = max(1, it->second.ki - gap);
        it->second.uc = it->second.ki + gap;
    }
    return fairnessConstraint;
}
