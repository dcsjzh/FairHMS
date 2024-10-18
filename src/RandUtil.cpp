#include "RandUtil.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

int Idum = 0;

double rand0()
{
    const int TAB_SIZE = 97;
    int j;
    static double y, v[TAB_SIZE];
    static bool iff = false;
    const double RAN_DIVISOR = double(RAND_MAX + 1UL);
    if (RAN_DIVISOR < 0)
    {
        cout << "RAN_DIVISOR " << RAN_DIVISOR << endl;
        exit(1);
    }
    if (Idum < 0 || !iff)
    {
        iff = true;
        srand(Idum);
        Idum = 1;
        for (j = 0; j < TAB_SIZE; j++)
        {
            rand();
        }
        for (j = 0; j < TAB_SIZE; j++)
        {
            v[j] = rand();
        }
        y = rand();
    }
    j = int(TAB_SIZE * (y / RAN_DIVISOR));
    y = v[j];
    v[j] = rand();
    return y / RAN_DIVISOR;
}

double randUnif(double lo, double hi)
{
    return rand0() * (hi - lo) + lo;
}

double randGauss()
{
    static bool iset = false;
    static double gset;
    if (!iset)
    {
        double v1, v2;
        double r = 2.0;
        while (r >= 1.0)
        {
            v1 = randUnif(-1, 1);
            v2 = randUnif(-1, 1);
            r = v1 * v1 + v2 * v2;
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

void RandUtil::get_random_direction(size_t dim, Point &rp)
{
    vector<double> coords;
    double len = 0;
    while (len == 0)
    {
        for (size_t i = 0; i < dim; i++)
        {
            double rand = randGauss();
            len += rand * rand;
            coords.push_back(rand);
        }
        len = sqrt(len);
    }
    for (size_t i = 0; i < dim; i++)
    {
        coords[i] = coords[i] / len;
    }
    Point p(dim, coords);
    rp = p;
    return;
}
