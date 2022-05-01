#include "sphere.h"

/*
*	construct I in Sphere (Step 2)
*/
point_set_t* construct_I(int dim, int k)
{
	double R = 2 * sqrt((double) dim);

	// we assume k >= dim
	if (k < dim)
	{
		printf("dimension cannot be larger than k.\n");
		exit(0);
	}

	// Case 2 of Sphere where I contains a single point
	// The best radius of I is 2 * \sqrt(d) theoretically, but in practice, other radius also works well in some cases, e.g., \sqrt(d)
	if (k < dim * dim)
	{
		point_t* p = alloc_point(dim);
		for (int i = 0; i < dim; i++)
			p->coord[i] = 2;

		point_set_t* I = alloc_point_set(1);
		I->points[0] = scale(R / calc_len(p), p);
		release_point(p);
		return I;
	}

	// Case 3 of Sphere
	// t is the number of hyperplanes in each family that are used to divide the facet of hypercube (t >= 0)
	int t = pow(k / dim * 1.0 / dim, 1.0 / (dim - 1)) - 1;
	int k_pi = dim * pow(t + 1, dim - 1.0); // size of I

	point_set_t* I = alloc_point_set(k_pi);
	int index = 0;

	// build the input
	Vvi input = build_input(t, dim);

	// obtain the output which is the catersian product of (dim-1) set of doubles
	Vvi output;
	Vi outputTemp;
	cart_product(output, outputTemp, input.begin(), input.end());

	// with the (dim - 1) coordinates in hands, we let the remaining one coordinate to be 1
	// there are dim number of choices
	point_t* p = alloc_point(dim);
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < output.size(); j++)
		{

			Vi temp_p = output.at(j);
			for (int l = 0; l < dim; l++)
			{
				if (l < i)
					p->coord[l] = temp_p.at(l);
				else if (l == i)
					p->coord[l] = 1; // the i-th dimension is set to be 1
				else
					p->coord[l] = temp_p.at(l - 1);
			}

			// put the resulting point in I, remenber to scale the point so that it lies on the sphere of radius 
			I->points[index++] = scale( R / calc_len(p), p);
		}
	}
	release_point(p);

	return I;
}

// The primitive operation: basic computation
// Compute an updated basic with by incoporating p to the old basic X
// X must have size dim + 1
void basisComputation(point_set_t*& X, point_t*& NN, point_t* p, point_t* x)
{
	if (!isViolated(x, NN, p))
		return;

	// If the basic X is empty
	if (X->numberOfPoints == 0)
	{
		X->points[0] = p;
		X->numberOfPoints++;

		NN = copy(p);
		return;
	}

	int d = p->dim;

	// If the basic X is not empty, incorporate p to obtain a new basic
	X->points[X->numberOfPoints] = p;
	X->numberOfPoints++;


	//outer loop
	while (true)
	{
		// Compute the projection on the affine space
		point_t* q_X = projectPointsOntoAffineSpace(X, x);

		int m = X->numberOfPoints - 1;
		if (m == 0)
			return;

		//d * (m+1) matrix
		vector<double> line(m + 1, 0);
		vector< vector<double> > A(d, line);
		point_set_t* X_tmp = alloc_point_set(m);

		double maxLambda = -1;
		int maxT = -1;
		//m bounding hyperplanes
		for (int t = 0; t < X->numberOfPoints; t++)
		{
			//initialize X_tmp
			for (int i = 0; i < m; i++)
			{
				if (i < t)
					X_tmp->points[i] = X->points[i];
				else
					X_tmp->points[i] = X->points[i + 1];
			}

			//set matrix
			for (int i = 0; i < d; i++) {
				for (int j = 0; j < m - 1; j++) {
					A[i][j] = X_tmp->points[j + 1]->coord[i] - X_tmp->points[0]->coord[i];
				}

				A[i][m - 1] = NN->coord[i] - q_X->coord[i];
				A[i][m] = NN->coord[i] - X_tmp->points[0]->coord[i];
			}

			//solve system of linear equations
			vector<double> sol(m);
			sol = gaussNtimesD(A);

			double lambda = sol[m - 1];

			//find the minimum non-negative lambda
			if ((isZero(lambda) || lambda > 0) && maxLambda < -0.5)
			{
				maxLambda = lambda;
				maxT = t;
			}
			else if (maxLambda > -0.5 && lambda > 0 && !isZero(lambda))
			{
				if (lambda < maxLambda)
				{
					maxLambda = lambda;
					maxT = t;
				}
			}

		}
		release_point_set(X_tmp, false);

		// X is indeed the P-basis
		if (maxLambda > 1 && !isZero(maxLambda - 1))
		{
			release_point(NN);
			NN = q_X;
			return;
		}

		//set new X
		for (int i = 0; i < m; i++)
		{
			if (i < maxT)
				X->points[i] = X->points[i];
			else
				X->points[i] = X->points[i + 1];
		}
		X->numberOfPoints--;

		//set new NN: NN = NN + (q_X - NN) * lambda 
		for (int i = 0; i < d; i++)
			NN->coord[i] = NN->coord[i] + (q_X->coord[i] - NN->coord[i]) * maxLambda;

		release_point(q_X);
	}
}

// The recursion call of the P-basic computation in "A combinatorial bound for linear programming and related problems"
void search_basis(point_set_t* P, point_set_t* T, point_t* NN_T, point_t* x, point_set_t*& X, point_t*& NN_X)
{
	int dim = P->points[0]->dim;
	int count;

	point_set_t* Q = alloc_point_set(P->numberOfPoints);

	Q->numberOfPoints = T->numberOfPoints;
	for (int i = 0; i < T->numberOfPoints; i++)
		Q->points[i] = T->points[i];
	count = Q->numberOfPoints;

	X = alloc_point_set(dim + 1);
	X->numberOfPoints = T->numberOfPoints;
	for (int i = 0; i < T->numberOfPoints; i++)
		X->points[i] = T->points[i];
	NN_X = copy(NN_T);
	
	// check the points in P one by one
	for (int i = 0; i < P->numberOfPoints; i++)
	{
		point_t* p = P->points[i];

		bool isPointInT = false;
		for (int j = 0; j < T->numberOfPoints; j++)
		{
			if (p == T->points[j])
			{
				isPointInT = true;
				break;
			}
		}
		if (isPointInT)
			continue;

		Q->points[count] = p;
		count++;
		Q->numberOfPoints++;

		// The violation test: check whether we can improve the distance by incorporate p 
		if (isViolated(x, NN_X, p))
		{
			basisComputation(X, NN_X, p, x);

			point_set_t* newX;
			point_t* newNN_X;

			// recursion
			search_basis(Q, X, NN_X, x, newX, newNN_X);

			X = newX;
			NN_X = newNN_X;
		}

	}

	T->numberOfPoints = dim + 1;
	release_point_set(T, false);
	release_point(NN_T);

	release_point_set(Q, false);
}

// The P-basic computation in "A combinatorial bound for linear programming and related problems"
point_set_t* search_basis(point_set_t* point_set_v, point_t* query)
{

	point_set_t* P = point_set_v;
	point_t* x = query;

	int dim = P->points[0]->dim;

	point_set_t* T = alloc_point_set(dim + 1);
	T->numberOfPoints = 1;
	T->points[0] = P->points[0];

	point_t* NN_T = copy(P->points[0]);

	point_set_t* X;
	point_t* NN_X;

	// recusion
	search_basis(P, T, NN_T, x, X, NN_X);

	release_point(NN_X);

	return X;
}

// The complete Sphere algorithm
point_set_t* sphereWSImpLP(point_set_t* point_set, int k)
{
	int dim = point_set->points[0]->dim;
	int N = point_set->numberOfPoints;

	point_set_t* result = alloc_point_set(k);

	// Step 1: Initialize to the boundary points
	double* b_value = new double[dim];
	for (int j = 0; j < dim; j++)							
	{
		b_value[j] = 0;
		for (int i = 0; i < point_set->numberOfPoints; i++)
		{
			if (point_set->points[i]->coord[j] > b_value[j])
			{
				b_value[j] = point_set->points[i]->coord[j];
				result->points[j] = point_set->points[i];
			}
		}
	}
	
	int count = dim;			
	if (k >= 2 * dim)		
	{

		// Step 2: Construct set I
		point_set_t* I;
		I = construct_I(dim, k - dim);

		// Step 3: For each point in I, search its P-basic
		for (int i = 0; i < I->numberOfPoints; i++)
		{
			point_set_t* basis = search_basis(point_set, I->points[i]);

			for (int j = 0; j < basis->numberOfPoints; j++)
			{
				//include the basic into the final solution
				bool isNew = true;
				for (int p = 0; p < count; p++)
				{
					if (basis->points[j]->id == result->points[p]->id)
					{
						isNew = false;
						break;
					}
				}

				if (isNew)
					result->points[count++] = basis->points[j];
			}

			basis->numberOfPoints = dim + 1;
			release_point_set(basis, false);
		}

		release_point_set(I, true);
	}

	if(count == k)
		return result;

	// Step 4: ImpGreedy
	double epsilon = 0.0000001;
	double worst, w;
	point_t* max;
	int* active = new int[point_set->numberOfPoints];

	point_t* lastRound_max;
	point_set_t* directions = alloc_point_set(N);
	double *rr = new double[N];

	// compute mrr_O(p) for each p in D
	max = NULL;
	worst = 0;
	for (int i = 0; i< N; i++)
	{
		active[i] = 1;

		directions->points[i] = alloc_point(dim);

		rr[i] = worstDirection(count, result, point_set->points[i], directions->points[i]->coord);

		w = rr[i];
		if (worst<w)
		{
			worst = w;
			max = point_set->points[i];
		}
		if (rr[i] <= 0.0 + epsilon)
			active[i] = 0;
	}

	result->points[count++] = max;
	lastRound_max = max;

	while (count < k)
	{
		// Find a point with maximum regret
		max = NULL;
		worst = 0;

		for (int i = 0; i < N; i++)
		{
			// upper bounding
			if (worst > rr[i])
			{
				directions->points[i]->id = -2; // means the worst utility vector is outdated
				continue;
			}
			// invariant checking fails or the worst utility vector is outdated
			else if (directions->points[i]->id == -2 || dot_prod(sub(point_set->points[i], lastRound_max), directions->points[i]) < rr[i])
			{
				// solve the exact LP
				rr[i] = worstDirection(count, result, point_set->points[i], directions->points[i]->coord);
				directions->points[i]->id = -1;
			}

			// update the current worst point
			w = rr[i];
			if (worst<w)
			{
				worst = w;
				max = point_set->points[i];
			}
			if (w <= 0.0 + epsilon)
				active[i] = 0;
		}

		// Add a point if regret > 0. Otherwise, stop.
		if (worst >= 0.0 + epsilon)
		{
			result->points[count++] = max;
			lastRound_max = max;
		}
		else
		{
			break;
		}
	}

	// fill in any remaining points with the first point
	for (int j = count; j < k; ++j)
		result->points[j] = result->points[0];


	delete[] rr;
	delete[] active;
	delete[] b_value;
	
	release_point_set(directions, true);
	return result;
}


void runSphere(vector<Point> dataP, int r, int k, vector<Point> &result, vector<Point> curSky, double &time)
{
	cout << "Running Sphere..." << endl;
	point_set_t *fatP4c, *result4c;
	fatP4c = RMSUtils::pointSetTransf(dataP);
	clock_t start = clock();
	result4c = sphereWSImpLP(fatP4c, r);
	clock_t end = clock();
	time = (end - start) / (double)(CLOCKS_PER_SEC / 1000);
	RMSUtils::pointSetTransf(result, result4c);
	release_point_set(fatP4c, true);
}