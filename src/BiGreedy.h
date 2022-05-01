#ifndef BIGREEDY_H
#define BIGREEDY_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <queue>
#include<map>
#include "Point.h"
#include"MyUtils.h"
#include "gurobi_c++.h"

void runBiGreedyAlg(unordered_map<int, vector<Point>> &categorizedDataP, 
                      unordered_map<int, fairinCate> &partitionMatroid,  
                      vector<Point> &dataP,                               
                      int cateId,                                        
                      double epsilon,                           
                      vector<Point> &allUtility,                          
                      vector<int> &result,
                      int k, int maxM, double &time);

void runBiGreedyPlusAlg(unordered_map<int, vector<Point>> &categorizedDataP, 
                      unordered_map<int, fairinCate> &partitionMatroid,          
                      vector<Point> &dataP,                               
                      int cateId,                                        
                      double gamma,                                      
                      double epsilon,                                    
                      vector<Point> &allUtility,                          
                      vector<int> &result,
                      int k, int maxM, double &time);

vector<int> constructCandidate(unordered_map<int, vector<Point>> &categorizedDataP,
                                      unordered_map<int, fairinCate> &partitionMatroid,   
                                      vector<int> &result,
                                      unordered_map<int, int> fairofR, 
                                      int k);

#endif