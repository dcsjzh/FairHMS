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
int elements_num, type_num;
int lower_range = 0, upper_range = 1, range_dis = upper_range - lower_range; 
map<string,int> type_map;  
map<int, int> type_total;

struct point 
{
    int ID;              
    double x, y;         
    double A_p0, dis_p0; 
    double A, B;         
    int type;            
}max_x, max_y;           
vector<struct point> element;
vector<struct point> convex_hull;                // points in a convex hull
struct envelope_point                            // points in an upper envelope
{
    double x, y, A, B;
    double Lx, Ly, Rx, Ry;
};
vector<struct envelope_point> envelope_points;   
struct cover                                     
{
    int ID;
    double L, R;
    friend bool operator < (struct cover const &a,struct cover const &b)
	{
        if(abs(a.L - b.L) < acc) 
        {
            if(abs(a.R - b.R) < acc)
                return a.ID < b.ID;
            return a.R < b.R;
        }
        return a.L < b.L;
	}
};
vector<set<struct cover>> cover_segments;   
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
vector<double> tau_set;     
struct cover_path           // structure of preserving the coverage point set for optimizing the DP algorithm
{
    bool visited = false;
    double cover = 0;
    vector<int> path;
    vector<int> types_count;
    bool operator<(const cover_path& a) const
    {
        if(abs(cover - a.cover) < acc)
        {
            int element_count_self = 0, element_count = 0;
            for(int i = 0; i < types_count.size(); i++)
                element_count_self += types_count[i];
            for(int i = 0; i < a.types_count.size(); i++)
                element_count += a.types_count[i];
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
bool compare_P0(point a, point b)
{
    if(abs(a.A_p0-b.A_p0) < acc)
        return a.dis_p0 > b.dis_p0;
    return a.A_p0 < b.A_p0 ;
}
bool compare_ID(point a, point b)
{
    return a.ID < b.ID;
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
int get_type(vector<string> data, vector<int> type_index)
{
    string type = "";
    for (int i = 0; i < type_index.size(); i++)
        type += data[type_index[i]];
    
    if(type_map.find(type) == type_map.end())
        type_map[type] = type_map.size();
    type_total[type_map[type]]++;
    return type_map[type];
}

void get_data(string file_path, vector<int> type_index)
{
    struct point ele;
    int i = 1;

	ifstream fin;
	fin.open(file_path, ios::in);
	string buff;
    getline(fin, buff);
	while (getline(fin, buff))
	{
        ele.ID = i++;
		vector<string> it = split(buff, '\t');
        double x = str2double(it[0]), y = str2double(it[1]);
        ele.type = get_type(it,type_index);
        ele.x = x;
        ele.y = y;
        ele.A = (x-y) / (range_dis);
        ele.B = y;
        ele.A_p0 = 0;
        ele.dis_p0 = 0;
        element.push_back(ele);
	}
	fin.close();
    elements_num = element.size();
    type_num = type_map.size();
}

/* prepare parameters required to calculate the upper envelope */
void new_elements()
{
    max_x = element[0];
    max_y = element[0];
    for(int i = 1; i < elements_num; i++)  
    {
        if(abs(max_x.x - element[i].x) < acc)
            max_x = max_x.y > element[i].y ? max_x : element[i];
        else
            max_x = max_x.x > element[i].x ? max_x : element[i];
        if(abs(max_y.y - element[i].y) < acc)
            max_y = max_y.x > element[i].x ? max_y : element[i];
        else
            max_y = max_y.y > element[i].y ? max_y : element[i];
    }
    for(int i = 0; i < elements_num; i++)
    {
        double x = max_x.x - element[i].x, y = max_x.y - element[i].y;
        if(abs(x) < acc)
        {
            element[i].A_p0 = 0;
            element[i].dis_p0 = 0;
            continue;
        }
        element[i].A_p0 = y / x;
        element[i].dis_p0 = x*x + y*y;
    }
    sort(element.begin(), element.end(), compare_P0);
}

void get_envelope()
{
    if(max_x.x == max_y.x && max_x.y == max_y.y) 
    {
        convex_hull.push_back(max_x);
        struct envelope_point P;
        P.x = max_x.x; P.y = max_x.y; P.A = max_x.A; P.B = max_x.B;
        P.Lx = 0; P.Ly = max_x.x; P.Rx = range_dis; P.Ry = max_x.y;
        envelope_points.push_back(P); 
        return;
    }

    convex_hull.push_back(max_x);
    convex_hull.push_back(element[0]);
    for(int i = 1; i < elements_num; i++)
    {
        int num = convex_hull.size();
        if(element[i].y <= convex_hull[num-1].y) continue;
        double x = convex_hull[num-2].x - convex_hull[num-1].x, y = convex_hull[num-2].y - convex_hull[num-1].y;
        double A = y/x, B = convex_hull[num-2].y - convex_hull[num-2].x * A;

        if(A * element[i].x + B > element[i].y)
            convex_hull.push_back(element[i]);
        else if(abs(A * element[i].x + B - element[i].y) < acc)
        {
            convex_hull.pop_back();
            convex_hull.push_back(element[i]);
        }
        else
        {
            convex_hull.pop_back();
            i--;
        }
        if(convex_hull.back().y == max_y.y) break;
    }
    vector<point> swap;
    for(int i = convex_hull.size() - 1; i >= 0; i--)
        swap.push_back(convex_hull[i]);
    convex_hull = swap;
    
    struct envelope_point P; 
    P.Lx = 0; P.Ly = convex_hull[0].B;
    int num = convex_hull.size();
    for(int i = 0; i < num-1; i++)
    {
        P.Rx = (convex_hull[i+1].B - convex_hull[i].B) / (convex_hull[i].A - convex_hull[i+1].A);
        P.Ry = convex_hull[i].A * P.Rx + convex_hull[i].B;
        P.x = convex_hull[i].x; P.y = convex_hull[i].y; P.A = convex_hull[i].A; P.B = convex_hull[i].B;
        envelope_points.push_back(P);
        P.Lx = P.Rx; P.Ly = P.Ry;
    }
    P.Rx = range_dis; P.Ry = convex_hull[num-1].x;
    P.x = convex_hull[num-1].x; P.y = convex_hull[num-1].y; P.A = convex_hull[num-1].A; P.B = convex_hull[num-1].B;
    envelope_points.push_back(P);
}

/* get the set composed of all possible \tau's */
void get_tau_set()
{
    vector<pair<double, double>> intersection;
    for(int i = 0; i < elements_num; i++)
    {
        intersection.push_back({0,element[i].x});
        intersection.push_back({1,element[i].y});
        for(int j = i + 1; j < elements_num; j++)
        {
            double x = (element[i].B - element[j].B)/(element[j].A - element[i].A);
            double y = element[i].A * x + element[i].B;
            if(x < 0 || x > range_dis) continue;
            intersection.push_back({x,y});
        }
    }
    sort(intersection.begin(), intersection.end(), compare_x);
    int j = 0;
    for(int i = 0; i < intersection.size(); i++)
    {
        while(intersection[i].first > envelope_points[j].Rx) j++;
        double y = envelope_points[j].A * intersection[i].first + envelope_points[j].B;
        tau_set.push_back(intersection[i].second/y);
    }
    sort(tau_set.begin(), tau_set.end());
}

void get_cover(double tau) 
{
    cover_segments.clear();
    for (int i = 0; i < type_num; i++)
        cover_segments.push_back({});    
    int envelope_num = envelope_points.size();
    for(int i = 0; i < elements_num; i++)
    {
        int j;
        struct cover Cov;

        for(j = 0; j<envelope_num && element[i].A*envelope_points[j].Lx+element[i].B < envelope_points[j].Ly*tau; j++);
        if(j == 0)
            Cov.L = 0;
        else if(j == envelope_num)
        {
            if(element[i].y >= envelope_points.back().Ry*tau)
            {
                Cov.L = (element[i].B - envelope_points.back().B*tau)/(envelope_points.back().A*tau - element[i].A);
                Cov.R = envelope_points.back().Rx;
                Cov.ID = element[i].ID;
                cover_segments[element[i].type].insert(Cov);
            }
            continue;
        }
        else
            Cov.L = (element[i].B - envelope_points[j-1].B*tau)/(envelope_points[j-1].A*tau - element[i].A);
        
        for(; j<envelope_num && element[i].A*envelope_points[j].Rx+element[i].B > envelope_points[j].Ry*tau; j++);
        if(j == envelope_num)
            Cov.R = envelope_points.back().Rx;
        else
            Cov.R = (element[i].B - envelope_points[j].B*tau)/(envelope_points[j].A*tau - element[i].A);
        
        Cov.ID = element[i].ID;
        cover_segments[element[i].type].insert(Cov);
    }
    DP_log.clear();
}

/* calculate the quota allocated to each group according to the proportions

vector<int> get_types_count(int num)
{
    vector<int> types_count;
    for(int i = 0; i < type_total.size(); i++)
    {
        int a = num * type_total[i] / (double)(elements_num) + 0.5;
        types_count.push_back(a);
    }
    return types_count;
}
*/

/* calculate the rightmost value of the coverage under the fair constraint by the DynProg algorithm */
pair<double, vector<int>> DP_cover_stack(vector<int> types_count_up, vector<int> types_count_down, int k)
{
    stack<struct cover_path> cover_stack;
    struct cover_path top_cover_path;
    top_cover_path.types_count = types_count_up;
    cover_stack.push(top_cover_path);
    pair<double, vector<int>> result = {0,{}};
    vector<int> zero_types;
    for (int i = 0; i < types_count_up.size(); i++)
        zero_types.push_back(0);
    DP_log[zero_types] = result;

    while(!cover_stack.empty())
    {
        if(cover_stack.top().visited)
        {
            top_cover_path = cover_stack.top();
            int total_k = 0;
            for(int i = 0; i < top_cover_path.types_count.size(); i++)
                total_k += max(top_cover_path.types_count[i], types_count_down[i]);
            if(total_k > k)
            {
                cover_stack.pop();
                continue;
            }
            for(int i = 0; i < types_count_up.size(); i++)
            {
                vector<int> old_path = top_cover_path.types_count;
                if(old_path[i] == 0) continue;
                old_path[i]--;
                pair<double, vector<int>> new_path = DP_log[old_path];

                struct cover Cov;
                vector<int> path = new_path.second;
                Cov.ID = 0;
                Cov.R = new_path.first;
                Cov.L = new_path.first;
                for(set<cover>::iterator it = cover_segments[i].begin(); it != cover_segments[i].end() && it->L - new_path.first < acc; it++)
                    if(Cov.R < it->R)
                        Cov = (*it);
                if(Cov.ID != 0)
                {
                    path.push_back(Cov.ID);
                    new_path = {Cov.R, path};
                }

                if(abs(top_cover_path.cover - new_path.first) < acc) 
                    top_cover_path.path = top_cover_path.path.size() < new_path.second.size() ? top_cover_path.path : new_path.second;
                else
                    if(top_cover_path.cover < new_path.first)
                    {
                        top_cover_path.cover = new_path.first;
                        top_cover_path.path = new_path.second;
                    }
                if(abs(top_cover_path.cover - range_dis) < acc) break;
            }
            pair<double, vector<int>> new_result = {top_cover_path.cover, top_cover_path.path};
            if(abs(new_result.first - range_dis) < acc)
                return new_result;
            DP_log[top_cover_path.types_count] = new_result;
            result = result.first > new_result.first ? result : new_result;
            cover_stack.pop();
        }
        else
        {
            cover_stack.top().visited = true;
            vector<int> top_types_count = cover_stack.top().types_count;
            for(int i = 0; i < top_types_count.size(); i++)
            {
                if(top_types_count[i] == 0) continue;
                vector<int> old_path = top_types_count;
                old_path[i]--;
                if(DP_log.find(old_path) != DP_log.end()) continue;
                struct cover_path new_cover_path;
                new_cover_path.types_count = old_path;
                cover_stack.push(new_cover_path);
            }
        }
        int asd = 0; 
        asd++;   
    }
    return result;
}

/* calculate the maximum tau by binary search */
pair<double, vector<int>> Binary_tau_DP(vector<int> types_count_up, vector<int> types_count_down, int k)
{
    int L = 0, R = tau_set.size();
    for(int M = (R-L)/2 + L; R-L > 1; M = (R-L)/2 + L)
    {
        get_cover(tau_set[M]);
        if(abs(DP_cover_stack(types_count_up, types_count_down, k).first - range_dis) < acc)
            L = M;
        else
            R = M;
    }
    get_cover(tau_set[L]);
    pair<double, vector<int>> result = {1 - tau_set[L], DP_cover_stack(types_count_up, types_count_down, k).second};

    return result;
}

/* randomly select points to satify the constraint */
pair<double, vector<int>> fill_N_acc(pair<double, vector<int>> result, vector<int> types_count)
{
    sort(element.begin(), element.end(), compare_ID);

    map<int, int> types_need;
    set<int> result_set;
    for(int i = 0; i < result.second.size(); i++)
    {
        types_need[element[result.second[i]-1].type]++;
        result_set.insert(result.second[i]);
    }
    for(int i = 0; i < types_count.size(); i++)
        types_need[i] = types_count[i] - types_need[i];
    for(int i = 0; i < elements_num; i++)
    {
        if(types_need[element[i].type] <= 0) continue;
        if(result_set.find(element[i].ID) == result_set.end())
        {
            result.second.push_back(element[i].ID);
            types_need[element[i].type]--;
        }
    }
    return result;
}

pair<double, vector<int>> total(string path, vector<int> type_index, vector<int> types_count_up, vector<int> types_count_down, int k, double &time)
{
    pair<double, vector<int>> result;
    get_data(path, type_index);
    new_elements();
    get_envelope();
    get_tau_set();

    clock_t start = clock();
    result = Binary_tau_DP(types_count_up, types_count_down, k);
    clock_t end = clock();
    time = end - start; 

    result = fill_N_acc(result, types_count_down);
    return result;
}

void runIntCov(string fileName, vector<int> &R, vector<int> cateIndex, vector<int> fairness_upper, vector<int> fairness_lower, int k, double &time, double &mhr)
{  
    acc = 0.000001;
    element.clear();
    type_map.clear();
    type_total.clear();
    convex_hull.clear();
    envelope_points.clear();
    cover_segments.clear();
    DP_log.clear();
    tau_set.clear();

    pair<double, vector<int>> result = total(fileName, cateIndex, fairness_upper, fairness_lower, k, time);
    mhr = 1.0 - result.first;

    R.clear();
    cout << result.first << endl;
    for (int i = 0; i < result.second.size(); i++)
        cout << result.second[i] << "\t";
    cout << endl;
    cout << "Time=" << time << endl;
    for(int id :result.second)
    {
        R.push_back(id-1);
    }

    for (int i = 0; i < result.second.size(); i++)
        cout << element[result.second[i]-1].type << "\t";
}
