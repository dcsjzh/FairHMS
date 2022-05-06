#include "Point.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include<iomanip>

Point::Point()
{
  
}

Point::Point(size_t dimension)
{
  this->dimension = dimension;
  this->id=-1;
  for(size_t i = 0; i < dimension; i++)
    coordinates.push_back(0);
}

Point::Point(size_t dimension, size_t id)
{
   this->dimension = dimension;
   this->id=id;
  for(size_t i = 0; i < dimension; i++)
    coordinates.push_back(0);
}

Point::Point(size_t dimension, vector<double>& coordinates)
{
  assert(dimension == coordinates.size());

  this->dimension = dimension;
  this->coordinates = coordinates;
}

Point::Point(size_t dimension, size_t id, vector<double>& coordinates)
{
  assert(dimension == coordinates.size());

  this->dimension = dimension;
  this->id = id;
  this->coordinates = coordinates;
}

Point::Point(size_t dimension, size_t id, vector<double>& coordinates, vector<int> &categorys)
{
  assert(dimension == coordinates.size());

  this->dimension = dimension;
  this->id = id;
  this->coordinates = coordinates;
  this->categorys = categorys;
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

void Point::set(size_t dimension, vector<double>& coordinates){
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

/*1 for this dominates other*/
bool Point::dominates(const Point& other_point) const
{
  assert(dimension == other_point.dimension);
  
  bool at_least_one = false;

  for(size_t i = 0; i < dimension; i++)
  {
    if(coordinates[i] < other_point.coordinates[i])
      return false;

    if(coordinates[i] >= other_point.coordinates[i])
      at_least_one = true;
  }

  return at_least_one;
}

void Point::dump(const char* prefix, const char* suffix) const
{
  cout << prefix;
  cout << "[ ";
  for(size_t i = 0; i < dimension-1; i++)
    cout << coordinates[i] << ", ";
  cout << coordinates[dimension-1];
  cout << " ]";
  cout << suffix;
  cout << endl;
}

Point Point::abs(const Point& other_point)
{
  Point p(other_point.dimension);
  for(size_t i = 0; i < other_point.dimension; i++)
    p.coordinates[i] = std::abs(other_point.coordinates[i]);

  return p;
}

Point Point::prod(const boost::numeric::ublas::matrix<double>& M, const Point& other_point)
{
  assert(M.size2() == other_point.dimension);
  Point p(M.size1());
  boost::numeric::ublas::vector<double> ublasp = to_ublas(other_point);
  ublasp = boost::numeric::ublas::prod(M,ublasp);
  p = from_ublas(ublasp);
  return p;
}

Point Point::from_ublas(const boost::numeric::ublas::vector<double>& ublasp)
{
  std::vector<double> coords;
  for(size_t i = 0; i < ublasp.size(); i++)
    coords.push_back(ublasp[i]);

  return Point(ublasp.size(), coords);
}

boost::numeric::ublas::vector<double> Point::to_ublas(const Point& p)
{
  boost::numeric::ublas::vector<double> ublasp(p.get_dimension());
  for(size_t i = 0; i < p.get_dimension(); i++)
    ublasp[i] = p.get_coordinate(i);

  return ublasp;
}

void Point::scale_to_length(double len)
{
  assert(len >= 0);

  double mylen = length();

  assert(len == 0 || mylen > 0);

  double factor = (mylen > 0) ? (len / mylen) : 0;

  for(size_t i = 0; i < dimension; i++)
    coordinates[i] *= factor;
  
}

void Point::print()
{
   std::cout<<"dim="<< dimension <<"  ";
    for(int i=0;i<dimension;i++)
    {
        std::cout <<setiosflags(ios::fixed)<<setprecision(6)<< coordinates[i] << " ";
    }
    std::cout<<endl;
}

vector<double> Point::getAllCoords()
{
  return this->coordinates;
}

double Point::dotP(vector<double> coord)
{
   double dotp = 0.0;
    for(size_t i = 0; i < dimension; i++)
      dotp += coordinates[i] * coord[i];

    return dotp;
}

Point& Point::operator=(const Point &p)
{
  //cout << "Point = "<<endl;
  if(this != &p)
  {
    dimension = p.dimension;
    id = p.id;
    coordinates = p.coordinates;
    categorys = p.categorys;
  }
  return *this;
}

int Point::get_category(int cateId){
  return categorys[cateId];
}