#include "IOUtil.h"
#include <assert.h>

void IOUtil::write_output_points(const char* fname,
                                 size_t dim,
                                 vector<Point>& dataP)
{
  FILE *fpt = fopen(fname,"w");

  if(!fpt)
  {
    cerr << "Cannot open file " << fname << " for writing \n";
    exit(1);    
  }

  if(dataP.size() == 0)
    goto done; 

  fprintf(fpt,"%ld\n",dim);

  for(size_t i = 0; i < dataP.size(); i++)
  {
    assert(dim == dataP[i].get_dimension());
    for(size_t j = 0; j < dim; j++)
    {
      fprintf(fpt,"%lf",dataP[i].get_coordinate(j));
      if(j < dim - 1)
        fprintf(fpt," ");
    }
    if(i < dataP.size() - 1)
      fprintf(fpt,"\n");
  }

done:
  fclose(fpt);
}

void IOUtil::read_input_points(const char* fname, 
		               size_t& dim, 
			       vector<Point>& dataP)
{

  FILE* fpt = fopen(fname, "r");

  if(!fpt)
  {
    cerr << "Cannot open file " << fname << " for reading \n";
    exit(1);    
  }

  fscanf(fpt,"%ld", &dim);

  double coord;
  size_t count=0;
  while(fscanf(fpt,"%lf",&coord) != EOF)
  {
    vector<double> coords;
    coords.push_back(coord);
    for(size_t i = 1; i < dim; i++)
    {
      fscanf(fpt,"%lf",&coord);
      coords.push_back(coord);
    }

    Point p(dim,count,coords);
    dataP.push_back(p);
    count++;
  }

  fclose(fpt);
}

void IOUtil::readUtils(const char* fname, size_t dim, vector<Point>& Utils, size_t size)
{
  FILE* fpt = fopen(fname, "r");

  if(!fpt)
  {
    cerr << "Cannot open file " << fname << " for reading \n";
    exit(1);    
  }
  for(size_t i=0;i<size;i++)
  {
    double coord;
    vector<double> coords;
    for(size_t i = 0; i < dim; i++)
    {
      fscanf(fpt,"%lf",&coord);
      coords.push_back(coord);
    }

    Point u(dim, i, coords);
    Utils.push_back(u);
  }
   fclose(fpt);
}

void IOUtil::readAllUtils(const char* fname, size_t dim, vector<Point>& Utils)
{
   FILE* fpt = fopen(fname, "r");

  if(!fpt)
  {
    cerr << "Cannot open file " << fname << " for reading \n";
    exit(1);    
  }
  double coord;
  int count=0;
  while(fscanf(fpt,"%lf",&coord) != EOF)
  {
    vector<double> coords;
    coords.push_back(coord);
    for(size_t i = 1; i < dim; i++)
    {
      fscanf(fpt,"%lf",&coord);
      coords.push_back(coord);
    }

    Point u(dim, count, coords);
    count++;
    Utils.push_back(u);
  }

  fclose(fpt);
}

void IOUtil::readRandomNumofUtils(const char* fname, size_t dim, vector<Point>& Utils, size_t size, int randomNumber)
{
  FILE* fpt = fopen(fname, "r");
  if(!fpt)
  {
    cerr << "Cannot open file " << fname << " for reading \n";
    exit(1);    
  }
  int count=0;
  //load qian d+1 utils
  for(size_t i=0;i<dim+1;i++)
  {
    double coord;
    vector<double> coords;
    for(size_t i = 0; i < dim; i++)
    {
      fscanf(fpt,"%lf",&coord);
      coords.push_back(coord);
    }

    Point u(dim,count, coords);
    Utils.push_back(u);
    count++;
  }
  //absorb random number of utils
  for(size_t i=0;i<randomNumber;i++)
  {
    double coord;
    vector<double> coords;
    for(size_t i = 0; i < dim; i++)
    {
      fscanf(fpt,"%lf",&coord);
      coords.push_back(coord);
    }
  }
  //load  fixed utils
  for(size_t i=0;i<size-dim-1;i++)
  {
    double coord;
    vector<double> coords;
    for(size_t i = 0; i < dim; i++)
    {
      fscanf(fpt,"%lf",&coord);
      coords.push_back(coord);
    }

    Point u(dim, count, coords);
    Utils.push_back(u);
    count++;
  }
   fclose(fpt);
}

void IOUtil::readRandomNumofUtilsNoD(const char* fname, size_t dim, vector<Point>& Utils, size_t size, int randomNumber)
{
  FILE* fpt = fopen(fname, "r");
  if(!fpt)
  {
    cerr << "Cannot open file " << fname << " for reading \n";
    exit(1);    
  }
  int count=0;
  //absorb random number of utils
  for(size_t i=0;i<randomNumber;i++)
  {
    double coord;
    vector<double> coords;
    for(size_t i = 0; i < dim; i++)
    {
      fscanf(fpt,"%lf",&coord);
      coords.push_back(coord);
    }
  }
  //load  fixed utils
  for(size_t i=0;i<size;i++)
  {
    double coord;
    vector<double> coords;
    for(size_t i = 0; i < dim; i++)
    {
      fscanf(fpt,"%lf",&coord);
      coords.push_back(coord);
    }

    Point u(dim, count, coords);
    Utils.push_back(u);
    count++;
  }
   fclose(fpt);
}

void IOUtil::readFixedSizePoints(const char* fname, size_t& dim, 
                                          vector<Point>& dataP, vector<int>& toInsert, size_t &size, vector<bool> &isDeleted)
{
  FILE* fpt = fopen(fname, "r");

  if(!fpt)
  {
    cerr << "Cannot open file " << fname << " for reading \n";
    exit(1);    
  }

  fscanf(fpt,"%ld", &dim);

//load init points
  double coord;
  size_t count=0;
  while(fscanf(fpt,"%lf",&coord) != EOF)
  {
    vector<double> coords;
    coords.push_back(coord);
    for(size_t i = 1; i < dim; i++)
    {
      fscanf(fpt,"%lf",&coord);
      coords.push_back(coord);
    }
    Point p(dim, count, coords);
    dataP.push_back(p);
    count++;
  }

  size = count / 2;
  for(int i=0;i<count;i++)
  {
    if(i <= (size - 1))
    {
        isDeleted.push_back(false);
    }
    else
    {
        toInsert.push_back(i);
        isDeleted.push_back(true);
    }
  }

  fclose(fpt);
}