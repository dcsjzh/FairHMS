#include "Point.h"
#include <cassert>
#include <cmath>
#include <iostream>
using namespace std;

Point::Point()
{

}

Point::Point(size_t dimension)
{
    this->dimension = dimension;
    for(size_t i = 0; i < dimension; i++)
    coordinates.push_back(0);
}

Point::Point(size_t dimension, vector<double>& coordinates)
{
    assert(dimension == coordinates.size());
    this->dimension = dimension;
    this->coordinates = coordinates;
}

Point Point::operator-(const Point& other_point) const
{
    assert(other_point.dimension == dimension);
    Point diff(dimension);
    for(size_t i = 0; i < dimension; i++)
    diff.coordinates[i] = coordinates[i] - other_point.coordinates[i];
    return diff;
}

Point Point::operator*(double factor) const
{
    Point mult(dimension);
    for(size_t i = 0; i < dimension; i++)
    mult.coordinates[i] = coordinates[i] * factor;
    return mult;
}

void Point::set(size_t dimension, vector<double>& coordinates)
{
    assert(dimension == coordinates.size());
    this->dimension = dimension;
    this->coordinates = coordinates;
}

double Point::distance_to(const Point& other_point) const
{
    Point diff = (*this) - other_point;
    return diff.length();
}

double Point::length() const
{
    double length = 0;
    for(size_t i = 0; i < dimension; i++)
    length += coordinates[i] * coordinates[i];
    return sqrt(length);
}

double Point::dotP(const Point& other_point) const
{
    assert(dimension == other_point.dimension);
    double dotp = 0.0;
    for(size_t i = 0; i < dimension; i++)
    dotp += coordinates[i] * other_point.coordinates[i];
    return dotp;
}

double Point::get_coordinate(size_t idx) const
{
    assert(idx < dimension);
    return coordinates[idx];
}

bool Point::dominates(const Point& other_point) const
{
    assert(dimension == other_point.dimension);
    bool at_least_one = false;
    for(size_t i = 0; i < dimension; i++)
    {
        if(coordinates[i] < other_point.coordinates[i])
        {
            return false;
        }
        if(coordinates[i] > other_point.coordinates[i])
        {
            at_least_one = true;
        }
    }
    return at_least_one;
}

void Point::dump(const char* prefix, const char* suffix) const
{
    cout << prefix;
    cout << "[ ";
    for(size_t i = 0; i < dimension-1; i++)
    {
        cout << coordinates[i] << ", ";
    }
    cout << coordinates[dimension-1];
    cout << " ]";
    cout << suffix;
    cout << endl;
}

Point Point::abs(const Point& other_point)
{
    Point p(other_point.dimension);
    for(size_t i = 0; i < other_point.dimension; i++)
    {
        p.coordinates[i] = std::abs(other_point.coordinates[i]);
    }
    return p;
}

void Point::scale_to_length(double len)
{
    assert(len >= 0);
    double mylen = length(); 
    assert(len == 0 || mylen > 0);
    double factor = (mylen > 0) ? (len / mylen) : 0;
    for(size_t i = 0; i < dimension; i++)
    {
        coordinates[i] *= factor;
        coordinates[i] *= coordinates[i];
    }  
}