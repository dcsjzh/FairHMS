/** 
 * This code has been downloaded from the website of Bojan Nikolic
 * see http://www.bnikolic.co.uk/blog/cpp-khachiyan-min-cov-ellipsoid.html
 * the code implements Khachiyan's algorithm to compute the Minimum volume
 * enclosing ellipsoid approximately. 
 *
 * Some parts were written inefficiently so have been changed.
 */

#include "MVE.h"
#include "RMSUtils.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "cholesky.h"

namespace ublas=boost::numeric::ublas;
template<class T>
bool InvertMatrix(const ublas::matrix<T> &input,
                    ublas::matrix<T> &inverse)
{
    using namespace boost::numeric::ublas;

    typedef permutation_matrix<std::size_t> pmatrix;
    matrix<T> A(input);
    pmatrix pm(A.size1());
    int res = lu_factorize(A,pm);
    if( res != 0 ) return false;
    inverse.assign(identity_matrix<T>(A.size1()));
    lu_substitute(A, pm, inverse);
    return true;
  }

void InvertLP(const ublas::matrix<double> &Lambdap,
                ublas::matrix<double> &LpInv)
{
    bool res=InvertMatrix(Lambdap, LpInv);
    if (not res)
    {
        throw std::runtime_error("Could not invert Matrix");
    }
}

void Lift(const ublas::matrix<double> &A,
            ublas::matrix<double> &Ap)
{
    Ap.resize(A.size1()+1,
              A.size2());
    ublas::matrix_range<ublas::matrix<double> >
      sub(Ap,
          ublas::range(0, A.size1()),
          ublas::range(0, A.size2()));
    sub.assign(A);
    ublas::row(Ap, Ap.size1()-1)=ublas::scalar_vector<double>(A.size2(),1.0);

}

void genDiag(const ublas::vector<double> &p,
               ublas::matrix<double> &res)
{
    res.assign(ublas::zero_matrix<double>(p.size(),
                                          p.size()));
    for(size_t i=0; i<p.size(); ++i)
    {
        res(i,i)=p(i);
    }
}

void KaLambda(const ublas::matrix<double> &Ap,
                const ublas::vector<double> &p,
                ublas::matrix<double> &Lambdap)
{
    //cout << "KaLambda: Ap.size1() = " << Ap.size1() << ", Ap.size2() = "
    //	 << Ap.size2() << ", p.size() = " << p.size() << endl;
    assert(p.size() == Ap.size2());

    /** 
        This code is very inefficient (It is trying to allocate a matrix
        of size m x m where m = # of points which is simply infeasible
        for large data sets

        It is basically trying to compute the sum of p_i q_i q_i^t for 
        i in 1 to m where m is the number of points and p is a m dimensional
        vector p = (p_1, ...., p_m)^t. Since each q_i is a n x 1 vector
        the return value of this is a n x n matrix.	
     **/
    /*
    ublas::matrix<double> dp(p.size(), p.size());
    genDiag(p, dp);

    dp=ublas::prod(dp, ublas::trans(Ap));
    Lambdap=ublas::prod(Ap,
                        dp);
    */

    /* Recall that Lambda(p) is a matrix of size n x n where it is simply
       defined as Lambda(p) = sum of p_i q_i q_i^t for i from 1 to m
       and q_i are the n x 1 vectors representing the (lifted points)
       The matrix Ap is (q_i .... q_m). It is a n x m matrix
    */
    size_t n = Ap.size1();
    ublas::matrix<double> sum(n,n);
    sum.assign(ublas::zero_matrix<double> (n, n));

    /* directly invoke the correct formula for matrix product */
    for(size_t i = 0; i < p.size(); i++)
    {
      for(size_t j = 0; j < n; j++)
        for(size_t k = 0; k < n; k++)
          sum(j,k) += p(i) * Ap(j, i) * Ap(k, i);
    }
    Lambdap = sum;
}


double KhachiyanIter(const ublas::matrix<double> &Ap,
                       ublas::vector<double> &p)
  {
    /// Dimensionality of the problem
    const size_t d=Ap.size1()-1;

    ublas::matrix<double> Lp;
    KaLambda(Ap, p, Lp);
    ublas::matrix<double> ILp(Lp.size1(), Lp.size2());
    InvertLP(Lp, ILp);


    /**
        * This code is very inefficient as it is trying to allocate
        * a matrix of size m x m where m is the number of data
        * points simply to get the diagonal elements of M !!

        * Basically it is trying to do the following:
	* M = ILp * Ap where ILp is a n x n matrix and
	* Ap is a n x m matrix (n is dimension, m is number
	* of points). The needed output is a set of m numbers
	* where the i-th number is basically <Ap_i^t, M_i>
	* i.e. the dot product of the transpose of the i-th
	* column of Ap and the i-th column of M. 

    ublas::matrix<double> M;
    M=ublas::prod(ILp, Ap);
    M=ublas::prod(ublas::trans(Ap), M);

    double maxval=0;
    size_t maxi=0;
    for(size_t i=0; i<M.size1(); ++i)
    {
        if (M(i,i) > maxval)
        {
            maxval=M(i,i);
            maxi=i;
        }
    }
    **/

    double maxval = 0;
    size_t maxi = 0;
    ublas::matrix<double> M = ublas::prod(ILp, Ap);
    ublas::matrix<double> Apt = ublas::trans(Ap);
    assert(M.size1() == Apt.size2());

    for(size_t i = 0; i < Ap.size2(); i++)
    {
      //evaluate the product c_i^t (ILp) c_i
      ublas::matrix_column<ublas::matrix<double> > mc(M, i);
      ublas::matrix_row<ublas::matrix<double> > mr(Apt,i);
     
      //evaluate the left product of row and column
      double prod = 0;
      for(size_t j = 0; j < mc.size(); j++)
        prod += mr(j) * mc(j);

      if(prod > maxval)
      {
        maxval = prod;
        maxi = i;
      }

    }

    const double step_size=(maxval -d - 1)/((d+1)*(maxval-1));
    ublas::vector<double> newp=p*(1-step_size);
    newp(maxi) += step_size;

    const double err= ublas::norm_2(newp-p);
    p=newp;
    return err;

}

void KaInvertDual(const ublas::matrix<double> &A,
                    const ublas::vector<double> &p,
                    ublas::matrix<double> &Q,
                    ublas::vector<double> &c
                    )
  {
    const size_t d=A.size1();
    /** This part of the code is not efficient. It is trying
        to allocate a matrix dp of size m x m where m is the
        number of points

	Basically it is trying to do the following: 
	Here A is a n x m matrix (where n is dimension
	and m is number of points)

	The output of this code is PN a m x n matrix where
	the row i is the transpose of the i-th column of A
	multiplied by p_i where p is a m-vector.
    **/
/*
    ublas::matrix<double> dp(p.size(), p.size());
    genDiag(p, dp);

    ublas::matrix<double> PN=ublas::prod(dp, ublas::trans(A));
*/
    //========= Begin replacement code ===============
    ublas::matrix<double> PN = ublas::trans(A);
    assert(p.size() == PN.size1());
    for(size_t i = 0; i < PN.size1(); i++)
      for(size_t j = 0; j < PN.size2(); j++)
        PN(i,j) *= p(i);
    //======== End replacement code ===================
    
    PN=ublas::prod(A, PN);

    ublas::vector<double> M2=ublas::prod(A, p);
    ublas::matrix<double> M3=ublas::outer_prod(M2, M2);

    ublas::matrix<double> invert(PN.size1(), PN.size2());
    InvertLP(PN - M3, invert);

    Q.assign( 1.0/d *invert);
    c=ublas::prod(A, p);
}


double KhachiyanAlgo(const ublas::matrix<double> &A,
                       double eps,
                       size_t maxiter,
                       ublas::matrix<double> &Q,
                       ublas::vector<double> &c)
{
    ublas::vector<double> p=ublas::scalar_vector<double>(A.size2(), 1.0)*(1.0/A.size2());

    ublas::matrix<double> Ap;
    Lift(A, Ap);

    double ceps=eps*2;
    for (size_t i=0;  i<maxiter && ceps>eps; ++i)
    {
        ceps=KhachiyanIter(Ap, p);
    }

    KaInvertDual(A, p, Q, c);

    return ceps;
}


void print_matrix(const ublas::matrix<double>& A)
{
   for(size_t i = 0; i < A.size1(); i++)
   {
     for(size_t j = 0; j < A.size2(); j++)
       cout << A(i,j) << "  ";
     cout << endl;
   }
}


void print_vector(const ublas::vector<double>& v)
{
  if(v.size() == 1)
  {
    cout << v(0) << endl;
    return;
  }

  cout << "---     ---" << endl;
  for(size_t i = 0; i < v.size(); i++)
  {
    cout << "||" << " " << v(i) << " " << "||" << endl;
  }

  cout << "---     ---" << endl;
}
void MVEUtil::GetNormalizedMVE(const vector<Point>& dataP,
                               float epsilon,
                               vector<Point>& normalizedP,
			       double& outer_rad,
			       double& inner_rad
                               )
{
    using namespace boost::numeric::ublas;

    if(dataP.size() == 0)
        return;
  
    size_t d = dataP[0].get_dimension();

    ublas::matrix<double> A(d,dataP.size());

    size_t j = 0;

    for(size_t i = 0; i < dataP.size(); i++)
    {
        for(size_t k=0; k < d; k++)
            A(k,j) = dataP[i].get_coordinate(k);
        ++j;
    }

    //try Khachiyans algorithm
    //cout << "Running Khachiyan algorithm " << endl;


    ublas::matrix<double> Q(d,d);
    ublas::vector<double> c(d);

    size_t maxiter = 1024 * dataP.size() * d;
    double ceps = KhachiyanAlgo(A,epsilon,maxiter,Q,c);
    if(ceps > epsilon)
    {
      throw std::runtime_error("Khachiyan failed ... ");
    }

    /** Only for debug    
    cout << "Khachiyan returned " << ceps << endl;
    cout << "Printing the matrix Q " << endl;
    print_matrix(Q);
    **/

    /**
        The equation of the ellipse that is returned
	by Khachiyan algorithm is (x - c)^T Q (x - c) <= 1
        Unfortunately all the points of the data set
	do not satisfy it (but they do satisfy it
	approximately). Now we can put (1+epsilon) * dim
	on the right hand side and the points will
	still satisfy it, but this may be too loose
	an approximation. Instead, we find the maximum
	of (x-c)^T Q (x-c) over all the data points,
	add 0.01, and take that as the RHS.
	The square root of this gives the radius of
	the outer sphere, when we apply a linear 
	transform to turn the ellipsoid into a ball
    **/
    double max_val = 0;
    for(size_t i = 0; i < dataP.size(); i++)
    {
      ublas::vector<double> ublasp = Point::to_ublas(dataP[i]);
      
      //compute p - c
      ublasp = ublasp - c;

      //compute Q (p - c)
      ublas::vector<double> prod1 = ublas::prod(Q,ublasp);

      //compute (p - c)^T
      ublas::matrix<double> ublasptr(1,d);
      for(size_t j = 0; j < d; j++)
        ublasptr(0,j) = ublasp[j];

      //compute (p-c)^T Q (p-c)
      ublas::vector<double> x = ublas::prod(ublasptr,prod1);

      if(x(0) > max_val) max_val = x(0);

      /** for debug
      cout << "Printing (p-c)^T Q (p -c ) " ;
      print_vector(x);
      cout << endl;
      **/
    }

    //the radius of the outer sphere
    outer_rad = sqrt(max_val + 0.01);


    //now compute the linear transformation
    ublas::matrix<double> L(d,d);
    //compute a Cholesky factorization
    int res = cholesky_decompose(Q,L);

    if(res != 0)
    {
      cout << "Cholesky decomposition failed " << res << endl;
      exit(1);
    }

    /** debug 
    cout << "Printing cholesky decomposition " << endl;
    print_matrix(L);
    **/

    //get the transpose of L
    ublas::matrix<double> LTr(d,d);
    LTr = ublas::trans(L);

    //ublas::matrix<double> chprod = ublas::prod(L, LTr);
    //print_matrix(chprod);

    normalizedP = dataP;

    for(size_t i = 0; i < dataP.size(); i++)
    {
      ublas::vector<double> ublasp = Point::to_ublas(normalizedP[i]);

      //compute the linear transformation
      ublasp = ublasp - c;
      ublasp = ublas::prod(LTr,ublasp);
      normalizedP[i] = Point::from_ublas(ublasp);
    }

    inner_rad = 1/((1+epsilon) * d);


    /* save this matrix into RandomUtil class for further use */
    /* the use requires the inverse of (matrix transpose) of Ltr */
    ublas::matrix<double> IL(d,d);
    InvertLP(L,IL);
    RMSUtils::dimension = d;
    RMSUtils::transformation_matrix = IL ;
    RMSUtils::center = c;

    return;
}

