#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <queue>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <stack>
using namespace std;

double acc = 0.00000001;   
int elements_num, group_num;
int lower_range = 0, upper_range = 1, range_dis = upper_range - lower_range; 
map<string,int> groups_ID;  

struct element 
{
    int idx;              
    double x, y;         
    double slope_rightest, dis_rightest; 
    double slope, intercept;         
    int group;            
}max_x, max_y;           
vector<struct element> ele_set;
vector<struct element> convex_Hull;                // points in a convex hull
struct envelope_ele                            // points in an upper envelope
{
    double x, y, slope, intercept;
    double Lx, Ly, Rx, Ry;
};
vector<struct envelope_ele> envelope;   
struct cover                                     
{
    int idx;
    double L, R;
    friend bool operator < (struct cover const &a,struct cover const &b)
	{
        if(abs(a.L - b.L) < acc) 
        {
            if(abs(a.R - b.R) < acc)
                return a.idx < b.idx;
            return a.R < b.R;
        }
        return a.L < b.L;
	}
};
vector<set<struct cover>> cover_Segments;   
struct HashFunc
{
	size_t operator() (const vector<int>& key) const {
		std::hash<int> hasher;
		size_t seed = 0;
		for (int i : key) {
			seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};
unordered_map<vector<int>,pair<double, vector<int>>, HashFunc> DP_log;  // store the intermediate states of a dynaimc program
vector<double> tau_candidate;     
struct cover_path           // structure of preserving the coverage point set for optimizing the DP algorithm
{
    bool visited = false;
    double cover = 0;
    vector<int> path;
    vector<int> groups_count;
    bool operator<(const cover_path& a) const
    {
        if(abs(cover - a.cover) < acc)
        {
            int element_count_self = 0, element_count = 0;
            for(int i = 0; i < groups_count.size(); i++)
                element_count_self += groups_count[i];
            for(int i = 0; i < a.groups_count.size(); i++)
                element_count += a.groups_count[i];
            return element_count_self > element_count;
        }
        return cover < a.cover; 
    }
};

double random_range(double x, double y)
{
    int max = 1000;
    double z, a = rand()%max, b = rand()%max;
    z = x + (a*max+b) / (max*max) * y;
    return z;
}
bool compare_rightest(element a, element b)
{
    if(abs(a.slope_rightest-b.slope_rightest) < acc)
        return a.dis_rightest > b.dis_rightest;
    return a.slope_rightest < b.slope_rightest ;
}
bool compare_idx(element a, element b)
{
    return a.idx < b.idx;
}
bool compare_x(pair<double,double> a,pair<double,double> b)
{
    return a.first < b.first;
}
static void _split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;

    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    _split(s, delim, elems);
    return elems;
}
double str2double(string num){
    double res;
    stringstream stream(num);
    stream>>res;
    return res;
}

/* transform the categrical attribute values into integers */
int get_group(vector<string> item, vector<int> group_col)
{
    string group = "";
    for (int i = 0; i < group_col.size(); i++)
        group += item[group_col[i]];
    
    if(groups_ID.find(group) == groups_ID.end())
        groups_ID[group] = groups_ID.size();
    return groups_ID[group];
}

void read_data_2D(string file_path, vector<int> group_col)
{
    struct element ele;
    int i = 1;

	ifstream fin;
	fin.open(file_path, ios::in);
	string buff;
    getline(fin, buff);
	while (getline(fin, buff))
	{
        ele.idx = i++;
		vector<string> it = split(buff, '\t');
        double x = str2double(it[0]), y = str2double(it[1]);
        ele.group = get_group(it,group_col);
        ele.x = x;
        ele.y = y;
        ele.slope = (x-y) / (range_dis);
        ele.intercept = y;
        ele.slope_rightest = 0;
        ele.dis_rightest = 0;
        ele_set.push_back(ele);
	}
	fin.close();
    elements_num = ele_set.size();
    group_num = groups_ID.size();
}

/* prepare parameters required to calculate the upper envelope */
void preprocess()
{
    max_x = ele_set[0];
    max_y = ele_set[0];
    for(int i = 1; i < elements_num; i++)  
    {
        if(abs(max_x.x - ele_set[i].x) < acc)
            max_x = max_x.y > ele_set[i].y ? max_x : ele_set[i];
        else
            max_x = max_x.x > ele_set[i].x ? max_x : ele_set[i];
        if(abs(max_y.y - ele_set[i].y) < acc)
            max_y = max_y.x > ele_set[i].x ? max_y : ele_set[i];
        else
            max_y = max_y.y > ele_set[i].y ? max_y : ele_set[i];
    }
    for(int i = 0; i < elements_num; i++)
    {
        double x = max_x.x - ele_set[i].x, y = max_x.y - ele_set[i].y;
        if(abs(x) < acc)
        {
            ele_set[i].slope_rightest = 0;
            ele_set[i].dis_rightest = 0;
            continue;
        }
        ele_set[i].slope_rightest = y / x;
        ele_set[i].dis_rightest = x*x + y*y;
    }
    sort(ele_set.begin(), ele_set.end(), compare_rightest);
}

void get_envelope()
{
    if(max_x.x == max_y.x && max_x.y == max_y.y) 
    {
        convex_Hull.push_back(max_x);
        struct envelope_ele ele;
        ele.x = max_x.x; ele.y = max_x.y; ele.slope = max_x.slope; ele.intercept = max_x.intercept;
        ele.Lx = 0; ele.Ly = max_x.x; ele.Rx = range_dis; ele.Ry = max_x.y;
        envelope.push_back(ele); 
        return;
    }

    convex_Hull.push_back(max_x);
    convex_Hull.push_back(ele_set[0]);
    for(int i = 1; i < elements_num; i++)
    {
        int num = convex_Hull.size();
        if(ele_set[i].y <= convex_Hull[num-1].y) continue;
        double x = convex_Hull[num-2].x - convex_Hull[num-1].x, y = convex_Hull[num-2].y - convex_Hull[num-1].y;
        double slope = y/x, intercept = convex_Hull[num-2].y - convex_Hull[num-2].x * slope;

        if(slope * ele_set[i].x + intercept > ele_set[i].y)
            convex_Hull.push_back(ele_set[i]);
        else if(abs(slope * ele_set[i].x + intercept - ele_set[i].y) < acc)
        {
            convex_Hull.pop_back();
            convex_Hull.push_back(ele_set[i]);
        }
        else
        {
            convex_Hull.pop_back();
            i--;
        }
        if(convex_Hull.back().y == max_y.y) break;
    }
    vector<element> reverse_CH;
    for(int i = convex_Hull.size() - 1; i >= 0; i--)
        reverse_CH.push_back(convex_Hull[i]);
    convex_Hull = reverse_CH;
    
    struct envelope_ele ele; 
    ele.Lx = 0; ele.Ly = convex_Hull[0].intercept;
    int num = convex_Hull.size();
    for(int i = 0; i < num-1; i++)
    {
        ele.Rx = (convex_Hull[i+1].intercept - convex_Hull[i].intercept) / (convex_Hull[i].slope - convex_Hull[i+1].slope);
        ele.Ry = convex_Hull[i].slope * ele.Rx + convex_Hull[i].intercept;
        ele.x = convex_Hull[i].x; ele.y = convex_Hull[i].y; ele.slope = convex_Hull[i].slope; ele.intercept = convex_Hull[i].intercept;
        envelope.push_back(ele);
        ele.Lx = ele.Rx; ele.Ly = ele.Ry;
    }
    ele.Rx = range_dis; ele.Ry = convex_Hull[num-1].x;
    ele.x = convex_Hull[num-1].x; ele.y = convex_Hull[num-1].y; ele.slope = convex_Hull[num-1].slope; ele.intercept = convex_Hull[num-1].intercept;
    envelope.push_back(ele);
}

/* get the set composed of all possible \tau's */
void get_tau_candidate()
{
    vector<pair<double, double>> intersections;
    for(int i = 0; i < elements_num; i++)
    {
        intersections.push_back({0,ele_set[i].x});
        intersections.push_back({1,ele_set[i].y});
        for(int j = i + 1; j < elements_num; j++)
        {
            double x = (ele_set[i].intercept - ele_set[j].intercept)/(ele_set[j].slope - ele_set[i].slope);
            double y = ele_set[i].slope * x + ele_set[i].intercept;
            if(x < 0 || x > range_dis) continue;
            intersections.push_back({x,y});
        }
    }
    sort(intersections.begin(), intersections.end(), compare_x);
    int j = 0;
    for(int i = 0; i < intersections.size(); i++)
    {
        while(intersections[i].first > envelope[j].Rx) j++;
        double y = envelope[j].slope * intersections[i].first + envelope[j].intercept;
        tau_candidate.push_back(intersections[i].second/y);
    }
    sort(tau_candidate.begin(), tau_candidate.end());
}

void get_cover(double tau) 
{
    cover_Segments.clear();
    for (int i = 0; i < group_num; i++)
        cover_Segments.push_back({});    
    int envelope_num = envelope.size();
    for(int i = 0; i < elements_num; i++)
    {
        int j;
        struct cover Cov;

        for(j = 0; j<envelope_num && ele_set[i].slope*envelope[j].Lx+ele_set[i].intercept < envelope[j].Ly*tau; j++);
        if(j == 0)
            Cov.L = 0;
        else if(j == envelope_num)
        {
            if(ele_set[i].y >= envelope.back().Ry*tau)
            {
                Cov.L = (ele_set[i].intercept - envelope.back().intercept*tau)/(envelope.back().slope*tau - ele_set[i].slope);
                Cov.R = envelope.back().Rx;
                Cov.idx = ele_set[i].idx;
                cover_Segments[ele_set[i].group].insert(Cov);
            }
            continue;
        }
        else
            Cov.L = (ele_set[i].intercept - envelope[j-1].intercept*tau)/(envelope[j-1].slope*tau - ele_set[i].slope);
        
        for(; j<envelope_num && ele_set[i].slope*envelope[j].Rx+ele_set[i].intercept > envelope[j].Ry*tau; j++);
        if(j == envelope_num)
            Cov.R = envelope.back().Rx;
        else
            Cov.R = (ele_set[i].intercept - envelope[j].intercept*tau)/(envelope[j].slope*tau - ele_set[i].slope);
        
        Cov.idx = ele_set[i].idx;
        cover_Segments[ele_set[i].group].insert(Cov);
    }
    DP_log.clear();
}

/* calculate the rightmost value of the coverage under the fair constraint by the DynProg algorithm */
pair<double, vector<int>> DP_cover_stack(vector<int> group_upper_bounds, vector<int> group_lower_bounds, int k)
{
    stack<struct cover_path> cover_stack;
    struct cover_path top_cover_path;
    top_cover_path.groups_count = group_upper_bounds;
    cover_stack.push(top_cover_path);
    pair<double, vector<int>> cover_result = {0,{}};
    vector<int> zero_groups;
    for (int i = 0; i < group_upper_bounds.size(); i++)
        zero_groups.push_back(0);
    DP_log[zero_groups] = cover_result;

    while(!cover_stack.empty())
    {
        if(cover_stack.top().visited)
        {
            top_cover_path = cover_stack.top();
            int total_k = 0;
            for(int i = 0; i < top_cover_path.groups_count.size(); i++)
                total_k += max(top_cover_path.groups_count[i], group_lower_bounds[i]);
            if(total_k > k)
            {
                cover_stack.pop();
                continue;
            }
            for(int i = 0; i < group_upper_bounds.size(); i++)
            {
                vector<int> tmp_cover_path = top_cover_path.groups_count;
                if(tmp_cover_path[i] == 0) continue;
                tmp_cover_path[i]--;
                pair<double, vector<int>> last_cover_path = DP_log[tmp_cover_path];

                struct cover Cov;
                vector<int> path = last_cover_path.second;
                Cov.idx = 0;
                Cov.R = last_cover_path.first;
                Cov.L = last_cover_path.first;
                for(set<cover>::iterator it = cover_Segments[i].begin(); it != cover_Segments[i].end() && it->L - last_cover_path.first < acc; it++)
                    if(Cov.R < it->R)
                        Cov = (*it);
                if(Cov.idx != 0)
                {
                    path.push_back(Cov.idx);
                    last_cover_path = {Cov.R, path};
                }

                if(abs(top_cover_path.cover - last_cover_path.first) < acc) 
                    top_cover_path.path = top_cover_path.path.size() < last_cover_path.second.size() ? top_cover_path.path : last_cover_path.second;
                else
                    if(top_cover_path.cover < last_cover_path.first)
                    {
                        top_cover_path.cover = last_cover_path.first;
                        top_cover_path.path = last_cover_path.second;
                    }
                if(abs(top_cover_path.cover - range_dis) < acc) break;
            }
            pair<double, vector<int>> new_result = {top_cover_path.cover, top_cover_path.path};
            if(abs(new_result.first - range_dis) < acc)
                return new_result;
            DP_log[top_cover_path.groups_count] = new_result;
            cover_result = cover_result.first > new_result.first ? cover_result : new_result;
            cover_stack.pop();
        }
        else
        {
            cover_stack.top().visited = true;
            vector<int> top_groups_count = cover_stack.top().groups_count;
            for(int i = 0; i < top_groups_count.size(); i++)
            {
                if(top_groups_count[i] == 0) continue;
                vector<int> tmp_cover_path = top_groups_count;
                tmp_cover_path[i]--;
                if(DP_log.find(tmp_cover_path) != DP_log.end()) continue;
                struct cover_path new_cover_path;
                new_cover_path.groups_count = tmp_cover_path;
                cover_stack.push(new_cover_path);
            }
        }
        int asd = 0; 
        asd++;   
    }
    return cover_result;
}

/* calculate the maximum tau by binary search */
pair<double, vector<int>> Bi_Search_tau(vector<int> group_upper_bounds, vector<int> group_lower_bounds, int k)
{
    int L = 0, R = tau_candidate.size();
    for(int M = (R-L)/2 + L; R-L > 1; M = (R-L)/2 + L)
    {
        get_cover(tau_candidate[M]);
        if(abs(DP_cover_stack(group_upper_bounds, group_lower_bounds, k).first - range_dis) < acc)
            L = M;
        else
            R = M;
    }
    get_cover(tau_candidate[L]);
    pair<double, vector<int>> tau_result = {tau_candidate[L], DP_cover_stack(group_upper_bounds, group_lower_bounds, k).second};

    return tau_result;
}

/* randomly select points to satify the constraint */
pair<double, vector<int>> post_Process(pair<double, vector<int>> tau_result, vector<int> groups_count)
{
    sort(ele_set.begin(), ele_set.end(), compare_idx);

    map<int, int> groups_absence;
    set<int> tmp_path;
    for(int i = 0; i < tau_result.second.size(); i++)
    {
        groups_absence[ele_set[tau_result.second[i]-1].group]++;
        tmp_path.insert(tau_result.second[i]);
    }
    for(int i = 0; i < groups_count.size(); i++)
        groups_absence[i] = groups_count[i] - groups_absence[i];
    for(int i = 0; i < elements_num; i++)
    {
        if(groups_absence[ele_set[i].group] <= 0) continue;
        if(tmp_path.find(ele_set[i].idx) == tmp_path.end())
        {
            tau_result.second.push_back(ele_set[i].idx);
            groups_absence[ele_set[i].group]--;
        }
    }
    return tau_result;
}

pair<double, vector<int>> get_result(string path, vector<int> group_col, vector<int> group_upper_bounds, vector<int> group_lower_bounds, int k, double &time)
{
    pair<double, vector<int>> tau_result;
    read_data_2D(path, group_col);
    preprocess();
    get_envelope();
    get_tau_candidate();

    clock_t start = clock();
    tau_result = Bi_Search_tau(group_upper_bounds, group_lower_bounds, k);
    clock_t end = clock();
    time = end - start; 

    tau_result = post_Process(tau_result, group_lower_bounds);
    return tau_result;
}

void runIntCov(string fileName, vector<int> &result, vector<int> group_col, vector<int> fairness_upper, vector<int> fairness_lower, int k, double &time, double &mhr)
{  
    acc = 0.000001;
    ele_set.clear();
    groups_ID.clear();
    convex_Hull.clear();
    envelope.clear();
    cover_Segments.clear();
    DP_log.clear();
    tau_candidate.clear();

    pair<double, vector<int>> tau_result = get_result(fileName, group_col, fairness_upper, fairness_lower, k, time);
    mhr = tau_result.first;

    result.clear();
    cout << tau_result.first << endl;
    for (int i = 0; i < tau_result.second.size(); i++)
        cout << tau_result.second[i] << "\t";
    cout << endl;
    cout << "Time=" << time << endl;
    for(int idx :tau_result.second)
    {
        result.push_back(idx-1);
    }

    for (int i = 0; i < tau_result.second.size(); i++)
        cout << ele_set[tau_result.second[i]-1].group << "\t";
}
