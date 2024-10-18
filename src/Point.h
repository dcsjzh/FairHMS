#ifndef __POINT_H__
#define __POINT_H__

#include <vector>

using namespace std;

class Point
{
    size_t dimension;
    vector<double> coordinates;

public:
    Point();
    Point(size_t dimension);
    Point(size_t dimension, vector<double> &coordinates);
    size_t get_dimension() const { return dimension; }
    vector<double> get_coordinates() const { return coordinates; }
    double get_coordinate(size_t idx) const;
    double dotP(const Point &other_point) const;
    bool dominates(const Point &other_point) const;
    void scale_to_length(double len);
    Point operator-(const Point &other_point) const;
    Point operator*(const double factor) const;
    void set(size_t dimension, vector<double> &coordinates);
    double distance_to(const Point &other_point) const;
    double length() const;
    void dump(const char *prefix, const char *suffix) const;
    static Point abs(const Point &other_point);
};

#endif
