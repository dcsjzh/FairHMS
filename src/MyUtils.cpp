#include "MyUtils.h"

vector<string> split(string src, char c)
{
    int i = 0;
    vector<string> ans;
    string tmp;
    while (i < src.length())
    {
        if (src[i] != c)
        {
            tmp += src[i++];
        }
        else
        {
            ans.push_back(tmp);
            tmp = "";
            ++i;
        }
    }
    ans.push_back(tmp);
    return ans;
}

void readDataPoint(string fileName, vector<Point> &dataP, int &dim, int &categoryDim,
                   vector<unordered_map<string, int>> &categoryToInt)
{
    ifstream fin(fileName);
    if (fin)
    {
        string line;
        // the first line represents the dimensionality
        getline(fin, line);
        dim = stoi(line);
        int maxNumofCategory = 5;
        vector<int> categoryCodes(maxNumofCategory, 0); 
        int count = 0;
        while (getline(fin, line))
        {
            vector<double> coords;
            vector<int> cates;
            vector<string> tmp = split(line, '\t');
            categoryDim = tmp.size() - dim;
            for (int j = 0; j < dim; ++j)
            {
                coords.push_back(stod(tmp[j]));
            }
            for (int j = 0; j < categoryDim; ++j)
            {
                if (categoryToInt[j].count(tmp[dim + j]) == 0)
                {
                    categoryToInt[j][tmp[dim + j]] = categoryCodes[j]++;
                }
                cates.push_back(categoryToInt[j][tmp[dim + j]]);
            }
            Point p(dim, count++, coords, cates);
            dataP.push_back(p);
        }
    }
    else
    {
        cout << "open \"" + fileName + "\" failed." << endl;
    }
}

/* read the utility functions from the file */
void readUtilityFunctions(string fileName, vector<Point> &utility, int m)
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
            Point u(coords.size(), count++, coords);
            utility.push_back(u);
            if (count == m)
            {
                break;
            }
        }
    }
}

/* generate groups */
void groupGen(unordered_map<int, vector<Point>> &categorizedDataP,
                   int categoryId, vector<Point> dataP)
{
    for (Point p : dataP)
    {
        int pCate = p.get_category(categoryId);
        if (categorizedDataP.count(pCate) == 0)
        {
            vector<Point> tmp;
            tmp.push_back(p);
            categorizedDataP[pCate] = tmp;
        }
        else
        {
            categorizedDataP[pCate].push_back(p);
        }
    }
}

/* write the results to the folder “result” */
void writeToFile(vector<Point> &dataP, int dim, vector<unordered_map<string, int>> &categoryToInt, unordered_map<int, fairinCate> &partitionMatroid,
                 string dataset, vector<int> &result, int k, int cateID, string algName, double time, double &MHR)
{
    // calculate the MHR of the final solution
    size_t ndir = RMSUtils::ndir_for_validation(dim);
    vector<Point> R;
    for (int pid : result)
    {
        R.push_back(dataP[pid]);
    }
    double MaxR = 1.0, AvgR, perc80;
    cout << "Computing HR ..." << endl;
    if (R.size() == 0)
    {
        MaxR = 0.0;
        MHR = - 1.0;
    }
    else
    {
        RMSUtils::Max_Avg_Regret(1.0, dim, ndir, dataP, R, MaxR, AvgR, perc80, 1);

        MHR = 1.0 - MaxR;
    }

    string outFileName = "../result/";
    outFileName += dataset;
    outFileName = outFileName.substr(0, outFileName.length() - 4);
    outFileName += "_" + to_string(k) + ".txt";
    ofstream fout(outFileName, ios::app);
    fout << algName << "categoryID=" << cateID << endl
         << "Fairness constraint: ";
    
    vector<string> cateNames;
    for (int j = 0; j < partitionMatroid.size(); ++j)
    {
        for (auto it = categoryToInt[cateID].begin(); it != categoryToInt[cateID].end(); ++it)
        {
            if (it->second == j)
            {
                cateNames.push_back(it->first);
                break;
            }
        }
    }
    
    for (int i = 0; i < partitionMatroid.size(); ++i)
    {
        fout << cateNames[i] << "=" << partitionMatroid[i].lc << "-" << partitionMatroid[i].ki << "-" << partitionMatroid[i].uc << " ";
    }
    fout << endl;
    vector<int> numofCates(partitionMatroid.size(), 0);
    for (int j = 0; j < partitionMatroid.size(); ++j)
    {
        int cateCode = j;
        string cateStr = cateNames[j];
        for (int pid : result)
        {
            if (dataP[pid].get_category(cateID) == cateCode)
            {
                fout << "idx=" << pid << "\t" << cateStr << "\t";
                for (int k = 0; k < dim; ++k)
                {
                    fout << to_string(dataP[pid].get_coordinate(k)) << "\t";
                }
                fout << endl;
                ++numofCates[cateCode];
            }
        }
    }
    
    for (int i = 0; i < partitionMatroid.size(); ++i)
    {
        fout << cateNames[i] << "=" << numofCates[i] << " ";
    }
    fout << endl;
    fout << "HR=" << 1.0 - MaxR << endl
         << "Time=" << time << endl
         << endl;
    fout.close();
}

/* write the results of the 2D algorithm - IntCov */
void writeToFile2D(vector<Point> &dataP, int dim, vector<unordered_map<string, int>> &categoryToInt, unordered_map<int, fairinCate> &partitionMatroid,
                   string dataset, vector<int> &result, int k, int cateID, string algName, double time, double mhr)
{
    string outFileName = "../result/";
    outFileName += dataset;
    outFileName = outFileName.substr(0, outFileName.length() - 4);
    outFileName += "_" + to_string(k) + ".txt";
    ofstream fout(outFileName, ios::app);
    fout << algName << "categoryID=" << cateID << endl
         << "Fairness constraint: ";
    
    vector<string> cateNames;
    for (int j = 0; j < partitionMatroid.size(); ++j)
    {
        for (auto it = categoryToInt[cateID].begin(); it != categoryToInt[cateID].end(); ++it)
        {
            if (it->second == j)
            {
                cateNames.push_back(it->first);
                break;
            }
        }
    }

    for (int i = 0; i < partitionMatroid.size(); ++i)
    {
        fout << cateNames[i] << "=" << partitionMatroid[i].lc << "-" << partitionMatroid[i].ki << "-" << partitionMatroid[i].uc << " ";
    }
    fout << endl;
    vector<int> numofCates(partitionMatroid.size(), 0);
    for (int j = 0; j < partitionMatroid.size(); ++j)
    {
        int cateCode = j;
        string cateStr = cateNames[j];
        for (int pid : result)
        {
            if (dataP[pid].get_category(cateID) == cateCode)
            {
                fout << "idx=" << pid << "\t" << cateStr << "\t";
                for (int k = 0; k < dim; ++k)
                {
                    fout << to_string(dataP[pid].get_coordinate(k)) << "\t";
                }
                fout << endl;
                ++numofCates[cateCode];
            }
        }
    }
 
    for (int i = 0; i < partitionMatroid.size(); ++i)
    {
        fout << cateNames[i] << "=" << numofCates[i] << " ";
    }
    fout << endl;
    fout << "HR=" << mhr << endl
         << "Time=" << time << endl
         << endl;
    fout.close();
}

/* generate the fairness constraints */
unordered_map<int, fairinCate> constructFairCons(unordered_map<int, vector<Point>> categorizedDataP, vector<Point> dataP, int k)
{
    unordered_map<int, fairinCate> partitionMatroid;
    
    vector<double> propotion(categorizedDataP.size());
    int total = 0, minCateId, minCateNum = k + 1;
    for (int i = 0; i < propotion.size(); ++i)
    {
        propotion[i] = categorizedDataP[i].size() / (dataP.size() * 1.0);
        partitionMatroid[i].ki = round(k * propotion[i]);
        total += partitionMatroid[i].ki;
        if (partitionMatroid[i].ki < minCateNum)
        {
            minCateNum = partitionMatroid[i].ki;
            minCateId = i;
        }
    }
    while (total < k)
    {
        ++partitionMatroid[minCateId].ki;
        ++total;
        minCateNum = partitionMatroid[minCateId].ki;
        //
        for (int i = 0; i < propotion.size(); ++i)
        {
            if (partitionMatroid[i].ki < minCateNum)
            {
                minCateNum = partitionMatroid[i].ki;
                minCateId = i;
            }
        }
    }
    
    int count = 0;
    for (auto it = partitionMatroid.begin(); it != partitionMatroid.end(); ++it)
    {
        if (it->second.ki == 0)
        {
            ++(it->second.ki);
            ++count;
        }
    }
    while (count != 0)
    {
        int maxCode, maxCount = 0;
        for (auto it = partitionMatroid.begin(); it != partitionMatroid.end(); ++it)
        {
            if (it->second.ki > maxCount)
            {
                maxCode = it->first;
                maxCount = it->second.ki;
            }
        }
        if (maxCount > 1)
        {
            --partitionMatroid[maxCode].ki;
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
    for (auto it = partitionMatroid.begin(); it != partitionMatroid.end(); ++it)
    {
        it->second.lc = max(1, it->second.ki - gap);
        it->second.uc = it->second.ki + gap;
    }
    return partitionMatroid;
}
