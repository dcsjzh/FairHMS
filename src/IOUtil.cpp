#include "IOUtil.h"
#include "Point.h"
#include <cassert>
#include <cmath>
#include <ctime>
using namespace std;

template <typename T>
T RandT(T _min, T _max)
{
    T temp;
    if (_min > _max)
    {
        temp = _max;
        _max = _min;
        _min = temp;
    }
    return rand() / (double)RAND_MAX * (_max - _min) + _min;
}

void IOUtil::read_input_points(const char *fname, size_t &dim, vector<Point> &dataP)
{
    FILE *fpt = fopen(fname, "r");
    if (!fpt)
    {
        cerr << "Cannot open file " << fname << " for reading \n";
        exit(1);
    }
    fscanf(fpt, "%ld", &dim);
    double coord;
    while (fscanf(fpt, "%lf", &coord) != EOF)
    {
        vector<double> coords;
        coords.push_back(coord);
        for (size_t i = 1; i < dim; i++)
        {
            fscanf(fpt, "%lf", &coord);
            coords.push_back(coord);
        }
        Point p(dim, coords);
        dataP.push_back(p);
    }
    fclose(fpt);
}

void IOUtil::read_input_functions(const char *fname, size_t &dim, vector<UtilityFunction> &FC)
{
    srand((unsigned)time(NULL));
    FILE *fpt = fopen(fname, "r");
    if (!fpt)
    {
        cerr << "Cannot open file " << fname << " for reading \n";
        exit(1);
    }
    fscanf(fpt, "%ld", &dim);
    double coord;
    while (fscanf(fpt, "%lf", &coord) != EOF)
    {
        vector<double> coords;
        coords.push_back(coord);
        for (size_t i = 1; i < dim; i++)
        {
            fscanf(fpt, "%lf", &coord);
            coords.push_back(coord);
        }
        Point p(dim, coords);
        UtilityFunction utilityfunction;
        utilityfunction.direction = p;
        utilityfunction.fmax = 0;
        FC.push_back(utilityfunction);
    }
    fclose(fpt);
}

void IOUtil::write_output_points(const char *fname, size_t dim, vector<Point> &dataP)
{
    FILE *fpt = fopen(fname, "w");
    if (!fpt)
    {
        cerr << "Cannot open file " << fname << " for writing \n";
        exit(1);
    }
    if (dataP.size() == 0)
    goto done;
    for (size_t i = 0; i < dataP.size(); i++)
    {
        for (size_t j = 0; j < dim; j++)
        {
            fprintf(fpt, "%lf", dataP[i].get_coordinate(j));
            if (j < dim - 1)
            fprintf(fpt, " ");
        }
        if (i < dataP.size() - 1)
        {
            fprintf(fpt, "\n");
        }
    }
    done:
    fclose(fpt);
}