#ifndef __POINT_H__
#define __POINT_H__

#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

class Point{
    size_t dimension;

    size_t id;

    vector<double> coordinates;

    vector<int> categorys;//

  public:
    Point();
    Point(size_t dimension);
    Point(size_t dimension, size_t id);
    Point(size_t dimension, vector<double>& coordinates);
    Point(size_t dimension, size_t id, vector<double>& coordinates);
    Point(size_t dimension, size_t id, vector<double>& coordinates, vector<int> &categorys);
    Point operator -(const Point& other_point) const;
    Point operator *(const double factor) const;
    void set(size_t dimension, vector<double>& coordinates);
    size_t get_dimension() const {return dimension;}
    size_t getId() const {return id;}
    double get_coordinate(size_t idx) const;
    double distance_to(const Point& other_point) const;
    double length() const;
    double dotP(const Point& other_point) const;
    bool dominates(const Point& other_point) const;
    void dump(const char* prefix, const char* suffix) const;
    void scale_to_length(double len); 
    static Point abs(const Point& other_point);    
    static Point prod(const boost::numeric::ublas::matrix<double>& M, const Point& other_point);
    static Point from_ublas(const boost::numeric::ublas::vector<double>& ublasp);
    static boost::numeric::ublas::vector<double> to_ublas(const Point& p); 
    void print();
    vector<double> getAllCoords();
    double dotP(vector<double> coord);

    Point& operator=(const Point &p);
    int get_category(int cateId);
    int getCategoryDim(){return categorys.size();}
};


#endif
