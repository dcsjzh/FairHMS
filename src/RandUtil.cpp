#include "RandUtil.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

/** 
    Global variables 
**/
int Idum = 0;

/**
    This is taken from the ANN library random number generator code
    by Sunil Arya and David Mount. They credit this function to
    William Press, Brian Flannery, Saul Teukolsky and William 
    Vetterling (see their book "Numerical Recipes in C").

    The function returns a uniform random number between 0 and 1
    using system supplied random() or rand(). The code does an
    additional randomizing shuffle to make it safer to use
**/
double rand0()
{
  const int TAB_SIZE = 97; //set to any large number

  int j;

  static double y, v[TAB_SIZE];

  static bool iff = false;

  const double RAN_DIVISOR = double(RAND_MAX + 1UL);

  if(RAN_DIVISOR < 0)
  {
    cout << "RAN_DIVISOR " << RAN_DIVISOR << endl;
    exit(1);
  }

  //Initialize on the first call, even if Idum is not set negative.
  //Determine MAXRAN, the next integer after the largest representable
  //value of type int. Assume this is a factor of 2 smaller than
  //the corresponding value of type unsigned int.
  if(Idum < 0 || !iff)
  {
    iff = true;
    srand(Idum); //reseed the generator

    Idum = 1;

    for(j = 0; j < TAB_SIZE; j++)
      rand(); //value intentionally ignored

    for(j = 0; j < TAB_SIZE; j++)
      v[j] = rand();
    y = rand();  //starting value
  }

  //if not initializing this is where we start. Use
  //previouslt saved random number y to get an index
  //between 1 and TAB_SIZE - 1. Then use corresponding
  //v[j] for both the next j and as output number

  j = int(TAB_SIZE * (y / RAN_DIVISOR));
  y = v[j];
  v[j] = rand(); //refill table entry

  return y / RAN_DIVISOR;

}

double randUnif(double lo, double hi)
{
  return rand0() * (hi - lo) + lo;
}


//return a gaussian random number with mean zero and variance 1
//Code taken from ann library 
double randGauss()
{
  static bool iset = false;
  static double gset;

  if(!iset)
  {
    double v1, v2;
    double r = 2.0;

    while(r >= 1.0)
    {
      //Pick two uniform numbers in the square [-1,1]^2 and see
      //if they are in the circle of radius 1. If not try again.
      v1 = randUnif(-1,1);
      v2 = randUnif(-1,1);
      r = v1*v1 + v2*v2;
    }
    double fac = sqrt(-2.0 * log(r) / r);
    gset = v1 * fac;
    iset = true;
    return v2 * fac;
  }
  else
  {

    iset = false;
    return gset;
  }
}

double RandUtil::randUnif(double lo, double hi)
{
  return randUnif(lo, hi);
}

void RandUtil::get_random_direction(size_t dim, Point& rp)
{
  vector<double> coords;

  //std::default_random_engine generator;
  //std::normal_distribution<double> distribution(0,1);
  double len = 0;
  while(len == 0)
  {
    for(size_t i = 0; i < dim; i++)
    {      
      double rand = randGauss();
     // double rand= distribution(generator);
      len += rand * rand;
      coords.push_back(rand);
    }

    len = sqrt(len);
  }
  for(size_t i = 0; i < dim; i++)
    coords[i] = coords[i] / len;
  Point p(dim,coords);
  rp = p;

/*
  cout << "[ ";
  for(size_t i = 0; i < dim; i++)
  {
    cout << p.get_coordinate(i);
    if(i < dim - 1)
      cout << ", ";
  }
  cout << " ]" << endl;
*/
  return; 
}

