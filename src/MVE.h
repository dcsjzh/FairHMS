#ifndef __MVE_H__
#define __MVE_H__

#include "Point.h"
#include <vector>

using namespace std;

class MVEUtil{
public:
  static void GetNormalizedMVE(const vector<Point>& dataP,
                        float epsilon,
                        vector<Point>& normalizedP,
			double& outer_rad,
			double& inner_rad
                       );

};

#endif
