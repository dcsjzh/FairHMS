#include "operation.h"

/*
 *	Check zero for floating points
 */
int isZero(double x)
{
	return x > -EQN_EPS && x < EQN_EPS;
}


/*
 *	Generate a random float value in [min_v, max_v].
 */
float rand_f( float min_v, float max_v)
{
	float rand_v;

	if( min_v > max_v)
		return 0;

	//srand( time(NULL));

	rand_v = float( rand( )) / RAND_MAX; //rand_v falls in [0,1] randomly.

	return ( max_v - min_v) * rand_v + min_v;
}

/*
 *	Calculate the Euclidean distance between two points
 */
DIST_TYPE calc_dist( point_t* point_v1, point_t* point_v2)
{
	int i, dim;
	DIST_TYPE diff;

	diff = 0;
	
	dim = point_v1->dim;
	for( i=0; i<dim; i++)
	{
		double t = point_v1->coord[ i];
		double t2 = point_v2->coord[ i];
		diff += ( DIST_TYPE)pow( point_v1->coord[ i] - point_v2->coord[ i], 2);
	}

	return ( DIST_TYPE)sqrt( diff);
}

/*
 *	Calculate the norm of a given vector
 */
DIST_TYPE calc_len(point_t* point_v)
{
	int i, dim;
	DIST_TYPE diff;

	diff = 0;
	
	dim = point_v->dim;
	for( i=0; i<dim; i++)
	{
		diff += point_v->coord[i] * point_v->coord[i];
	}

	return ( DIST_TYPE)sqrt( diff);

}


/*
 *  Create a new point give the target point
 *  (useful when some point need to be freed)
 */
point_t* copy(point_t * point_v2)
{
	if(point_v2 == NULL)
		return NULL;

	point_t* point_v1 = alloc_point(point_v2->dim);

	for(int i = 0; i< point_v2->dim; i++)
		point_v1->coord[i] = point_v2->coord[i];

	return point_v1;
}

/*
 *	Calculate the dot product between two points
 */
double dot_prod(point_t* point_v1, point_t* point_v2)
{
	int dim = point_v1->dim;
	double result = 0;

	int i;
	for(i = 0; i < dim; i++)
	{
		result += point_v1->coord[i]*point_v2->coord[i];
	}
	return result;
}

/*
 *	Calculate the dot product between two points
 */
double dot_prod(point_t* point_v1, double* v)
{
	int dim = point_v1->dim;
	double result = 0;

	int i;
	for(i = 0; i < dim; i++)
	{
		result += point_v1->coord[i]*v[i];
	}
	return result;
}

/*
 *	Calculate the subtraction between two points.
 */
point_t* sub(point_t* point_v1, point_t* point_v2)
{
	point_t* result = alloc_point(point_v1->dim);

	int i;
	for(i = 0; i < point_v1->dim; i++)
	{
		result->coord[i] = point_v1->coord[i] - point_v2->coord[i];
	}
	return result;
}

/*
 *	Calculate the addition between two points.
 */
point_t* add(point_t* point_v1, point_t* point_v2)
{
	point_t* result = alloc_point(point_v1->dim);

	int i;
	for(i = 0; i < point_v1->dim; i++)
	{
		result->coord[i] = point_v1->coord[i] + point_v2->coord[i];
	}
	return result;
}

/*
 *	Scale the given point
 */
point_t* scale(double c, point_t* point_v)
{
	point_t* result = alloc_point(point_v->dim);

	int i;
	for(i = 0; i < point_v->dim; i++)
	{
		result->coord[i] = point_v->coord[i] * c;
	}

	return result;
}

/*
 *  Violation test for basic computation
 *	Check whether point e violates the hyperplane
 *  which has the vector (normal_q - normal_p) as normal vector and pass through normal_p
 */
bool  isViolated(point_t* normal_q, point_t* normal_p, point_t* e)
{
	if(isZero(calc_dist(normal_q,normal_p)))
	{
		return true;
	}

	bool result;
	point_t* temp_normal = sub(normal_q, normal_p);
	point_t* temp = sub(e, normal_p);
	
	//use the dot product to determin the above/below relation
	if(dot_prod(temp_normal, temp) > 0 && !isZero(dot_prod(temp_normal, temp)))
		result = true;
	else
		result = false;

	release_point(temp_normal);
	release_point(temp);

	return result;
}

/*
 *  Find the point in p farthest in the direction v
 */
point_t* maxPoint(point_set_t* p, double *v)
{
	int N = p->numberOfPoints;

	int i, maxIndex = 0;
	double max = 0.0;
	
	for(i = 0; i < N; ++i)
		if (dot_prod(p->points[i], v) > max)
		{
			maxIndex = i;
			max = dot_prod(p->points[i], v);
		}

	return p->points[maxIndex];
}

/*
 *  Gauss Elimination
 */
vector<double> gaussNtimesD(vector< vector<double> > A) {
	int n = A.size();
	int d = A[0].size() - 1;


	for (int i = 0; i< d; i++) {
		// Search for maximum in this column
		double maxEl = abs(A[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<n; k++) {
			if (abs(A[k][i]) > maxEl) {
				maxEl = abs(A[k][i]);
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k<d + 1; k++) {
			double tmp = A[maxRow][k];
			A[maxRow][k] = A[i][k];
			A[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k<n; k++) {
			double c = -A[k][i] / A[i][i];
			for (int j = i; j<d + 1; j++) {
				if (i == j) {
					A[k][j] = 0;
				}
				else {
					A[k][j] += c * A[i][j];
				}
			}
		}
	}

	//check validity
	int count = 0;
	for (int i = 0; i < n; i++)
	{
		bool allZero = true;
		for (int j = 0; j < d + 1; j++)
		{
			if (!isZero(A[i][j]))
			{
				allZero = false;
				break;
			}
		}

		if (!allZero)
			count++;
	}
	if (count != d)
	{
		//printf("No Solution.\n");
		vector<double> x(d);
		x[d - 1] = -1;
		return x;
	}

	// Solve equation Ax=b for an upper triangular matrix A
	vector<double> x(d);
	
	for (int i = d - 1; i >= 0; i--) {
		x[i] = A[i][d] / A[i][i];
		for (int k = i - 1; k >= 0; k--) {
			A[k][d] -= A[k][i] * x[i];
		}
	}
	
	return x;
}

/*
 *  Project point on to an affine space
 */
point_t* projectPointsOntoAffineSpace(point_set_t* space, point_t* p)
{
	if (space->numberOfPoints < 1)
	{
		printf("ERROR in pt to affine\n");
		exit(0);
	}
	if (space->numberOfPoints == 1)
		return space->points[0];

	int dim = space->points[0]->dim;
	int n = space->numberOfPoints - 1;
	point_set_t* dirVecs = alloc_point_set(n);
	for (int i = 0; i < n; i++)
		dirVecs->points[i] = sub(space->points[i + 1], space->points[0]);

	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			// p_i = p_i - sum_j( dot(p_i, p_j) / |p_j|^2 * p_j)
			double c_j = dot_prod(dirVecs->points[i], dirVecs->points[j]) / dot_prod(dirVecs->points[j], dirVecs->points[j]);

			for (int k = 0; k < dim; k++)
				dirVecs->points[i]->coord[k] = dirVecs->points[i]->coord[k] - c_j * dirVecs->points[j]->coord[k];
		}

	}
	
	for (int i = 0; i < dirVecs->numberOfPoints; i++)
	{
		double norm = calc_len(dirVecs->points[i]);
		for (int j = 0; j < dim; j++)
			dirVecs->points[i]->coord[j] = dirVecs->points[i]->coord[j] / norm;
	}

	point_t* tmp = sub(p, space->points[0]);

	point_t* coord = alloc_point(n);
	for (int i = 0; i < n; i++)
		coord->coord[i] = dot_prod(tmp, dirVecs->points[i]);

	for (int i = 0; i < dim; i++)
	{
		tmp->coord[i] = 0;
		for (int j = 0; j < n; j++)
			tmp->coord[i] += coord->coord[j] * dirVecs->points[j]->coord[i];
	}

	point_t* proj = add(tmp, space->points[0]);

	release_point(tmp);
	release_point(coord);
	release_point_set(dirVecs, true);

	return proj;
}

/*
*	build the input for catersian product (in the form of vector of vector of double)
*  t is the number of hyperplanes in each family that are used to divide the facet of hypercube (t >= 0)
*  e.g. when dim = 3, t = 0, means no division and each facet has 1 hypercube
*		when dim = 3, t = 1, means each 2-d facet is divided by 2 hyperplan, one for each dimension, resulting in 4 hypercubes on each facet
*/
Vvi build_input(int t, int dim) {

	// the distance between concecutive parallel hyperplanes
	double dist_bet = 1.0 / (t + 1);

	Vvi vvi;
	Vi vi;
	for (int i = 0; i < t + 1; i++)
	{
		//make the points at the center of each hypercube
		vi.push_back(i * dist_bet + dist_bet / 2);
	}

	// the query points are put on the front facet of hypercubes
	// and therefore, there exists one dimension for each point with value 1.
	// we need another (dim - 1) number of points to define a point in dim-dimensional space
	// and each of this dimension are formed by the same set of value
	for (int i = 0; i < dim - 1; i++)
		vvi.push_back(vi);
	return vvi;
}

/*
*  codes adapted from stack over flow.
*  use recursion to compue the catersian product in the format of vector of vector of double.
*/
void cart_product(
	Vvi& rvvi,  // final result
	Vi&  rvi,   // current result 
	Vvi::const_iterator me, // current input
	Vvi::const_iterator end) // final input
{
	if (me == end) {
		// terminal condition of the recursion. We no longer have
		// any input vectors to manipulate. Add the current result (rvi)
		// to the total set of results (rvvvi).
		rvvi.push_back(rvi);
		return;
	}

	// need an easy name for my vector-of-ints
	const Vi& mevi = *me;
	for (Vi::const_iterator it = mevi.begin();
	it != mevi.end();
		it++) {
		// final rvi will look like "a, b, c, ME, d, e, f"
		// At the moment, rvi already has "a, b, c"
		rvi.push_back(*it);  // add ME
		cart_product(rvvi, rvi, me + 1, end); //add "d, e, f"
		rvi.pop_back(); // clean ME off for next round
	}
}

// read points from the input files
point_set_t* read_points(char* input)
{
	FILE* c_fp;

	char filename[MAX_FILENAME_LENG];
	sprintf(filename, "%s", input);

	if ((c_fp = fopen(filename, "r")) == NULL)
	{
		fprintf(stderr, "Cannot open the data file %s.\n", filename);
		exit(0);
	}

	int number_of_points, dim;
	fscanf(c_fp, "%i%i", &number_of_points, &dim);

	point_set_t* point_set = alloc_point_set(number_of_points);

	// read points line by line
	for (int i = 0; i < number_of_points; i++)
	{
		point_t* p = alloc_point(dim, i);
		for (int j = 0; j < dim; j++)
		{
			fscanf(c_fp, "%lf", &p->coord[j]);
		}
		point_set->points[i] = p;
	}

	fclose(c_fp);
	return point_set;
}

// The dominance in skyline computation
int dominates(point_t* p1, point_t* p2)
{
	int i;

	for (i = 0; i < p1->dim; ++i)
		if (p1->coord[i] < p2->coord[i])
			return 0;

	return 1;
}

// Compute the set of skyline points
point_set_t* skyline_point(point_set_t *p)
{
	int i, j, dominated, index = 0, m;
	point_t* pt;

	int* sl = new int[p->numberOfPoints];

	for (i = 0; i < p->numberOfPoints; ++i)
	{
		dominated = 0;
		pt = p->points[i];

		// check if pt is dominated by the skyline so far   
		for (j = 0; j < index && !dominated; ++j)
			if (dominates(p->points[sl[j]], pt))
				dominated = 1;

		if (!dominated)
		{
			// eliminate any points in current skyline that it dominates
			m = index;
			index = 0;
			for (j = 0; j < m; ++j)
				if (!dominates(pt, p->points[sl[j]]))
					sl[index++] = sl[j];

			// add this point as well
			sl[index++] = i;
		}
	}

	point_set_t* skyline = alloc_point_set(index);
	for (int i = 0; i < index; i++)
		skyline->points[i] = p->points[sl[i]];

	delete[] sl;
	return skyline;
}

// compute the orthope set
void insertOrth(double* &points, int &count, point_t* v)
{
	int dim = v->dim;
	int orthNum = pow(2.0, dim) - 1;

	for(int i = 0; i< dim; i++)
	{
		points[count * dim + i] = v->coord[i];
	}
	count++;



	int startEnumPow;
	for(int i=0; i< orthNum - 1; i++){
		startEnumPow = i+1; //from 1 to n-1
		for(int j=0; j<dim; j++){
			points[count * dim + j] = v->coord[j] * (startEnumPow % 2);
			startEnumPow /= 2;
			//printf("%lf ", startProjData[i][j]);
		}
		count++;
	}


}

