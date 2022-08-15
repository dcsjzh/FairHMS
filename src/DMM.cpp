#include "DMM.h"

point_set_t *discretize(int gamma, int dim)
{
	int size = pow(gamma, (double)dim - 1);
	point_set_t *F = alloc_point_set(size);
	double alpha = PI / 2 / gamma;

	int *theta = new int[dim - 1];
	for (int i = 0; i < dim - 1; i++)
	{
		theta[i] = 1;
	}

	for (int i = 0; i < size; i++)
	{
		//for (int i = 0; i < dim - 1; i++)
		//{
		//	printf("%d ", theta[i]);
		//}
		//printf("\n");
		double r = 1;
		point_t *pt = alloc_point(dim);

		for (int j = dim - 1; j >= 1; j--)
		{
			double angle = (theta[j - 1] - 0.5) * alpha;
			pt->coord[j] = r * cos(angle);
			double tmp = sin(angle);
			r = r * tmp;
		}
		pt->coord[0] = r;

		F->points[i] = pt;

		for (int j = 0; j < dim - 1; j++)
		{
			if (theta[j] < gamma)
			{
				theta[j]++;
				break;
			}
			else
				theta[j] = 1;
		}
	}

	delete[] theta;

	return F;
}

point_set_t *MRST_oracle(double **M, int n, int m, double eps, point_set_t *point_set)
{
	int **M_tmp = new int *[n];
	for (int i = 0; i < n; i++)
	{
		M_tmp[i] = new int[m];
	}

	int *count = new int[n];
	//printf("----eps: %lf----\n", eps);
	for (int i = 0; i < n; i++)
	{
		count[i] = 0;
		for (int j = 0; j < m; j++)
		{
			M_tmp[i][j] = (M[i][j] < eps || isZero(M[i][j] - eps)) ? 1 : 0;
			count[i] += M_tmp[i][j];
		}
	}

	int totalCount = m;
	vector<point_t *> selected;
	while (totalCount > 0)
	{
		int maxCount = 0;
		int maxI;
		for (int i = 0; i < n; i++)
		{
			if (maxCount < count[i])
			{
				maxI = i;
				maxCount = count[i];
			}
		}

		totalCount -= maxCount;
		selected.push_back(point_set->points[maxI]);

		for (int i = 0; i < m; i++)
		{
			if (M_tmp[maxI][i] == 1)
			{
				for (int j = 0; j < n; j++)
				{
					if (M_tmp[j][i] == 1)
					{
						M_tmp[j][i] = 0;
						count[j]--;
					}
				}
			}
		}
	}

	point_set_t *S = alloc_point_set(selected.size());

	for (int i = 0; i < selected.size(); i++)
	{
		S->points[i] = selected[i];
	}

	delete[] count;
	for (int i = 0; i < n; i++)
	{
		delete[] M_tmp[i];
	}
	delete[] M_tmp;

	return S;
}

point_set_t *DMM(point_set_t *point_set, int k)
{

	int dim = point_set->points[0]->dim;
	int gamma = 5;

	point_set_t *F = discretize(gamma, dim);

	int n = point_set->numberOfPoints;
	int m = F->numberOfPoints;

	double **M = new double *[n];
	for (int i = 0; i < n; i++)
	{
		M[i] = new double[m];
	}

	for (int j = 0; j < m; j++)
	{
		double max = 0;
		for (int i = 0; i < point_set->numberOfPoints; i++)
		{
			double temp = dot_prod(F->points[j], point_set->points[i]);
			if (temp > max)
				max = temp;
		}

		for (int i = 0; i < n; i++)
		{

			M[i][j] = 1 - dot_prod(F->points[j], point_set->points[i]) / max;
		}
	}

	std::set<double> sortedV_set;
	std::set<double>::iterator it;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			double value = M[i][j];
			sortedV_set.insert(M[i][j]);
		}
	}
	vector<double> sortedV;
	for (it = sortedV_set.begin(); it != sortedV_set.end(); ++it)
	{
		sortedV.push_back((*it));
	}

	int low = 0, high = sortedV.size();
	while (low < high)
	{

		int mid = (low + high) / 2;

		point_set_t *tmp = MRST_oracle(M, n, m, sortedV[mid], point_set);

		if (tmp->numberOfPoints <= k)
		{
			high = mid;
		}
		else
		{
			low = mid + 1;
		}
		release_point_set(tmp, false);
	}
	point_set_t *result = MRST_oracle(M, n, m, sortedV[low], point_set);

	for (int i = 0; i < n; i++)
	{
		delete[] M[i];
	}
	delete[] M;

	return result;
}

point_set_t *DMM_Greedy(point_set_t *point_set, int k)
{
	int dim = point_set->points[0]->dim;
	int gamma = 10;

	point_set_t *F = discretize(gamma, dim);

	int n = point_set->numberOfPoints;
	int m = F->numberOfPoints;

	double **M = new double *[n];
	for (int i = 0; i < n; i++)
	{
		M[i] = new double[m];
	}

	for (int j = 0; j < m; j++)
	{
		double max = 0;
		for (int i = 0; i < point_set->numberOfPoints; i++)
		{
			double temp = dot_prod(F->points[j], point_set->points[i]);
			if (temp > max)
				max = temp;
		}

		for (int i = 0; i < n; i++)
		{

			M[i][j] = 1 - dot_prod(F->points[j], point_set->points[i]) / max;
		}
	}

	if (k < dim)
	{
		printf("k < dim. \n");
		exit(0);
	}

	point_set_t *result = alloc_point_set(k);

	int *activeI = new int[n];
	for (int i = 0; i < n; i++)
		activeI[i] = 1;

	double *b_value = new double[dim];
	int *b_index = new int[dim];

	for (int j = 0; j < dim; j++)
	{
		b_value[j] = 0;
		for (int i = 0; i < n; i++)
		{
			if (point_set->points[i]->coord[j] > b_value[j])
			{
				b_value[j] = point_set->points[i]->coord[j];
				b_index[j] = i;
			}
		}
	}

	for (int i = 0; i < dim; i++)
	{
		result->points[i] = point_set->points[b_index[i]];
		activeI[b_index[i]] = 0;
	}
	delete[] b_value;
	delete[] b_index;

	int count = dim;

	while (count < k)
	{
		double max = -1;
		int maxJ = -1;
		for (int j = 0; j < m; j++)
		{
			int minI = -1;
			double min = INF;
			for (int i = 0; i < n; i++)
			{
				if (activeI[i])
					continue;

				if (M[i][j] < min)
				{
					min = M[i][j];
					minI = i;
				}
			}

			if (max < min)
			{
				max = min;
				maxJ = j;
			}
		}

		double selectMin = INF;
		int selectI = 0;
		for (int i = 0; i < n; i++)
		{
			if (!activeI[i])
				continue;

			if (selectMin > M[i][maxJ])
			{
				selectMin = M[i][maxJ];
				selectI = i;
			}
		}

		result->points[count++] = point_set->points[selectI];
		activeI[selectI] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		delete[] M[i];
	}
	delete[] M;
	delete[] activeI;

	return result;
}

void runDMMGreedy(vector<Point> dataP, int r, int k, vector<Point> &result, vector<Point> curSky, double &time)
{
	cout << "Running DMMGreedy..." << endl;
	point_set_t *fatP4c, *result4c;
	fatP4c = RMSUtils::pointSetTransf(dataP);
	clock_t start = clock();
	result4c = DMM_Greedy(fatP4c, r);
	clock_t end = clock();
	time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);
	RMSUtils::pointSetTransf(result, result4c);
	release_point_set(fatP4c, true);
	//   release_point_set(result4c, true);
}

void runDMMRRMS(vector<Point> dataP, int r, int k, vector<Point> &result, vector<Point> curSky, double &time)
{
	cout << "Running DMMRRMS ..." << endl;
	point_set_t *fatP4c, *result4c;
	fatP4c = RMSUtils::pointSetTransf(dataP);
	clock_t start = clock();
	result4c = DMM(fatP4c, r);
	clock_t end = clock();
	time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);
	RMSUtils::pointSetTransf(result, result4c);
	release_point_set(fatP4c, true);
	//   release_point_set(result4c, true);
}