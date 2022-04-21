# FairHMS
A C++ implementation of FairHMS in our paper "Happiness Maximizing Sets under Group Fairness Constraints". All baseline algorithms we evaluate in our experiments are also available online. Please refer to https://users.cs.duke.edu/~ssintos/kRMS_SEA/ for the implementation of eps-Kernel, Greedy, and HittingSet; and https://www.cse.ust.hk/~raywong/code/sigmod18-sphere.zip for the implementation of DMM-Greedy, DMM-RRMS, GeoGreedy, and Sphere.

# Usage
1. compile: cd src; make;
2. run: cd executable; ./run.out dataset NumberofCategories k is2D (e.g. "./run.out anti_6_100_skyline.txt 1 10 0" for high dimensionality, "./run.out anti_2_100_skyline.txt 1 5 1" for 2D);
3. The results are in the folder "results".
