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

void runBiGreedyAlg(unordered_map<int, vector<Point>> &groupedDataP, 
                      unordered_map<int, fairofGroup> &fairnessConstraint,  
                      vector<Point> &dataP,                               
                      int groupID,                                        
                      double epsilon,                           
                      vector<Point> &allUtility,                          
                      vector<int> &result,
                      int k, int maxM, double &time);

void runBiGreedyPlusAlg(unordered_map<int, vector<Point>> &groupedDataP, 
                      unordered_map<int, fairofGroup> &fairnessConstraint,          
                      vector<Point> &dataP,                               
                      int groupID,                                        
                      double gamma,                                      
                      double epsilon,                                    
                      vector<Point> &allUtility,                          
                      vector<int> &result,
                      int k, int maxM, double &time);


void runBiGreedyAlgWithDelta(unordered_map<int, vector<Point>> &groupedDataP, 
                            unordered_map<int, fairofGroup> &fairnessConstraint,  
                            vector<Point> &dataP,                               
                            int groupID,                                        
                            double epsilon,                           
                            vector<Point> &allUtility,                          
                            vector<int> &result,
                            int k, int maxM, double &time, double deltaC);


void runBiGreedyPlusAlgWithDelta(unordered_map<int, vector<Point>> &groupedDataP, 
                                unordered_map<int, fairofGroup> &fairnessConstraint,   
                                vector<Point> &dataP,                               
                                int groupID,                                         
                                double lambda,                                      
                                double epsilon,                                     
                                vector<Point> &utiFunClass,                          
                                vector<int> &result,
                                int k, int maxM, double &time, double deltaC);

vector<int> constructCandidate(unordered_map<int, vector<Point>> &groupedDataP,
                                      unordered_map<int, fairofGroup> &fairnessConstraint,   
                                      vector<int> &result,
                                      unordered_map<int, int> fairofResult, 
                                      int k);

#endif