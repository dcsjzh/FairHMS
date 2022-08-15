#include "lp.h"


//#define DEBUG_LP

// Takes an array of points s (of size N) and  a point pt and returns
// the direction in which pt has the worst regret ratio (in array v)
// as well as this regret ratio itself. 

/* We solve the following LP with col variables v[0], v[1], ..., v[D], x 
   and row variables q_1 ... q_{K-1}, r1, r_2
   
   Max x
   s.t. - (pt.a[0]-s[0].a[0])v[0] - (pt.a[1]-s[0].a[1])v[1] - ... - (pt.a[D-1]-s[0].a[D-1])v[D-1] + x = q_1
        - (pt.a[0]-s[1].a[0])v[0] - (pt.a[1]-s[1].a[1])v[1] - ... - (pt.a[D-1]-s[1].a[D-1])v[D-1] + x = q_2
        ...
        - (pt.a[0]-s[K-1].a[0])v[0] - (pt.a[1]-s[K-1].a[1])v[1] - ... - (pt.a[D-1]-s[K-1].a[D]-1)v[D-1] + x= q_K
           pt.a[0]v[0] + pt.a[1]v[1]  ....  + pt.a[D-1]v[D-1] = r1
          -pt.a[0]v[0] - pt.a[1]v[1]  ....  - pt.a[D-1]v[D-1] = r2
   variables have the following bounds
       0 <= v[j] < infty
       0 <= x < infty
       -infty < q_i  <=0
       -infty < r1 <= 1
       -infty < r2 <= -1
*/

double worstDirection(point_set_t *s, point_t* pt, double* &v)
{
	int K = s->numberOfPoints;
	int D = pt->dim;

    int* ia = new int[1+(K+5)*(D+5)];  //TODO: delete
	int* ja = new int[1+(K+5)*(D+5)];  //TODO: delete
    double* ar = new double[1+(K+5)*(D+5)];   //TODO: delete
    int i, j;   
    double epsilon=0.0000000000001; 


    glp_prob *lp;
    lp=glp_create_prob(); 
    glp_set_prob_name(lp, "max_regret_ratio");
    glp_set_obj_dir(lp, GLP_MAX);
    
 
    glp_add_rows(lp, K+2);  // add K+2 rows: q_1...q_k and r_1 and r_2
    // Add rows q_1 ... q_K 
    for (i=1; i<=K; i++) {
        char buf[10]; 
        sprintf(buf, "q%d", i);
        glp_set_row_name(lp, i, buf); 
        glp_set_row_bnds(lp, i, GLP_UP, 0.0, 0.0); // -infity < qi<=0
    }
       // Add rows r_1 and r_2
    glp_set_row_name(lp, K+1, "r1");  
    glp_set_row_bnds(lp, K+1, GLP_UP, 0.0, 1.0+epsilon); // r1 <= 1
    glp_set_row_name(lp, K+2, "r2");  
    glp_set_row_bnds(lp, K+2, GLP_UP, 0.0, -1.0+epsilon); // r2 <=-1
    
    
    glp_add_cols(lp, D+1);    // add D+1 columns: v[1] ... v[D] and x
    // Add col v[1] ... v[D]
    for (i=1; i<=D; i++) {
        char buf[10]; 
        sprintf(buf, "v%d", i);
        
        glp_set_col_name(lp, i, buf); 
        glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0); // 0 <= v[i] < infty
        glp_set_obj_coef(lp, i, 0.0);  // objective: maximize x ONLY
    }
    
    // Add col x
    glp_set_col_name(lp, D+1, "x");
    glp_set_col_bnds(lp, D+1, GLP_FR, 0.0, 0.0); // -infty <= x <= infty
    glp_set_obj_coef(lp, D+1, 1.0);  // objective: maximize x

    /*
            v[1]           v[2]   ...        v[D]    x
      q1 -(p-p_i)[0] -(p-p_i)[1] ... -(p-p_i)[D-1]   1
      q2  ......  same as q1  .................
      ...
      qK  ......  same as q1 .................
      r1   p[0]           p[1]    ...        p[D-1]  0
      r2  -p[0]         -p[1]    ...       -p[D-1]   0
    */
    
    int counter=1; 
    // set value on row q1 ... qk
    for (i=1; i<=K; i++) {
        
        #ifdef DEBUG_LP
        fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, i, D+1, 1.0); //DEBUG
        #endif
        
        ia[counter]=i; ja[counter]=D+1; ar[counter++]=1.0; // a["qi", "x"] =1
        for (j=1; j<=D; j++) {
            #ifdef DEBUG_LP
            fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, i, j, -(pt.a[j-1]-s[i-1].a[j-1])); //DEBUG
            #endif
            ia[counter]=i; ja[counter]=j; 
			ar[counter++] = -(pt->coord[j-1]-s->points[i-1]->coord[j-1]); // a["qi", "v[j]"] = -(pt-s[i-1])[j-1]
       }    
    }

    // set value on row r1 and r2
    #ifdef DEBUG_LP
    fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+1, D+1, 0); //DEBUG
    fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+2, D+1, 0); //DEBUG
    #endif
    ia[counter]=K+1; ja[counter]=D+1; ar[counter++]=0.0; // a["r1", "x"]=0 
    ia[counter]=K+2; ja[counter]=D+1; ar[counter++]=0.0; // a["r2", "x"]=0

    for (i=1; i<=D; i++){
        #ifdef DEBUG_LP
        fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+1, i, pt.a[i-1]); //DEBUG
        fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+2, i, -pt.a[i-1]); //DEBUG
        #endif
		ia[counter]=K+1; ja[counter]=i; ar[counter++]=pt->coord[i-1];   // e.g. a["r1", "v[1]"]=pt[0];   
		ia[counter]=K+2; ja[counter]=i; ar[counter++]=-pt->coord[i-1];   // e.g. a["r2", "v[1]"]=-pt[0];   
    }
    
    

    // loading data  
    glp_load_matrix(lp, counter-1, ia, ja, ar);    
    
    
    // Use this to print out the LP if you're debugging
    #ifdef DEBUG_LP
    glp_write_lp(lp, NULL, "testlp.lp");  // DEBUG
    #endif
    
    
    // running simplex
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev=GLP_MSG_OFF; // turn off all message by glp_simplex 
    glp_simplex(lp, &parm);    
    
    // write solution 
    #ifdef DEBUG_LP
    glp_print_sol(lp, "testlp.sol");  // DEBUG
    #endif
    
    // get values
    double regret_ratio=glp_get_obj_val(lp); 
    for (i=0; i<D; i++)
        v[i]=glp_get_col_prim(lp, i+1); // v[0] is at column 1, v[1] at col 2 ...
    
    
    glp_delete_prob(lp); // clean up
    delete []ia;
	delete []ja;
	delete []ar;
	
    return regret_ratio; 
}    
double worstDirection(int index, point_set_t *s, point_t* pt, double* &v)
{
	int K = index;
	int D = pt->dim;

    int* ia = new int[1+(K+5)*(D+5)];  //TODO: delete
	int* ja = new int[1+(K+5)*(D+5)];  //TODO: delete
    double* ar = new double[1+(K+5)*(D+5)];   //TODO: delete
    int i, j;   
    double epsilon=0.0000000000001; 


    glp_prob *lp;
    lp=glp_create_prob(); 
    glp_set_prob_name(lp, "max_regret_ratio");
    glp_set_obj_dir(lp, GLP_MAX);
    
 
    glp_add_rows(lp, K+2);  // add K+2 rows: q_1...q_k and r_1 and r_2
    // Add rows q_1 ... q_K 
    for (i=1; i<=K; i++) {
        char buf[10]; 
        sprintf(buf, "q%d", i);
        glp_set_row_name(lp, i, buf); 
        glp_set_row_bnds(lp, i, GLP_UP, 0.0, 0.0); // -infity < qi<=0
    }
       // Add rows r_1 and r_2
    glp_set_row_name(lp, K+1, "r1");  
    glp_set_row_bnds(lp, K+1, GLP_UP, 0.0, 1.0+epsilon); // r1 <= 1
    glp_set_row_name(lp, K+2, "r2");  
    glp_set_row_bnds(lp, K+2, GLP_UP, 0.0, -1.0+epsilon); // r2 <=-1
    
    
    glp_add_cols(lp, D+1);    // add D+1 columns: v[1] ... v[D] and x
    // Add col v[1] ... v[D]
    for (i=1; i<=D; i++) {
        char buf[10]; 
        sprintf(buf, "v%d", i);
        
        glp_set_col_name(lp, i, buf); 
        glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0); // 0 <= v[i] < infty
        glp_set_obj_coef(lp, i, 0.0);  // objective: maximize x ONLY
    }
    
    // Add col x
    glp_set_col_name(lp, D+1, "x");
    glp_set_col_bnds(lp, D+1, GLP_FR, 0.0, 0.0); // -infty <= x <= infty
    glp_set_obj_coef(lp, D+1, 1.0);  // objective: maximize x

    /*
            v[1]           v[2]   ...        v[D]    x
      q1 -(p-p_i)[0] -(p-p_i)[1] ... -(p-p_i)[D-1]   1
      q2  ......  same as q1  .................
      ...
      qK  ......  same as q1 .................
      r1   p[0]           p[1]    ...        p[D-1]  0
      r2  -p[0]         -p[1]    ...       -p[D-1]   0
    */
    
    int counter=1; 
    // set value on row q1 ... qk
    for (i=1; i<=K; i++) {
        
        #ifdef DEBUG_LP
        fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, i, D+1, 1.0); //DEBUG
        #endif
        
        ia[counter]=i; ja[counter]=D+1; ar[counter++]=1.0; // a["qi", "x"] =1
        for (j=1; j<=D; j++) {
            #ifdef DEBUG_LP
            fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, i, j, -(pt.a[j-1]-s[i-1].a[j-1])); //DEBUG
            #endif
            ia[counter]=i; ja[counter]=j; 
			ar[counter++] = -(pt->coord[j-1]-s->points[i-1]->coord[j-1]); // a["qi", "v[j]"] = -(pt-s[i-1])[j-1]
       }    
    }

    // set value on row r1 and r2
    #ifdef DEBUG_LP
    fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+1, D+1, 0); //DEBUG
    fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+2, D+1, 0); //DEBUG
    #endif
    ia[counter]=K+1; ja[counter]=D+1; ar[counter++]=0.0; // a["r1", "x"]=0 
    ia[counter]=K+2; ja[counter]=D+1; ar[counter++]=0.0; // a["r2", "x"]=0

    for (i=1; i<=D; i++){
        #ifdef DEBUG_LP
        fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+1, i, pt.a[i-1]); //DEBUG
        fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+2, i, -pt.a[i-1]); //DEBUG
        #endif
		ia[counter]=K+1; ja[counter]=i; ar[counter++]=pt->coord[i-1];   // e.g. a["r1", "v[1]"]=pt[0];   
		ia[counter]=K+2; ja[counter]=i; ar[counter++]=-pt->coord[i-1];   // e.g. a["r2", "v[1]"]=-pt[0];   
    }
    
    

    // loading data  
    glp_load_matrix(lp, counter-1, ia, ja, ar);    
    
    
    // Use this to print out the LP if you're debugging
    #ifdef DEBUG_LP
    glp_write_lp(lp, NULL, "testlp.lp");  // DEBUG
    #endif
    
    
    // running simplex
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev=GLP_MSG_OFF; // turn off all message by glp_simplex 
    glp_simplex(lp, &parm);    
    
    // write solution 
    #ifdef DEBUG_LP
    glp_print_sol(lp, "testlp.sol");  // DEBUG
    #endif
    
    // get values
    double regret_ratio=glp_get_obj_val(lp); 
    for (i=0; i<D; i++)
        v[i]=glp_get_col_prim(lp, i+1); // v[0] is at column 1, v[1] at col 2 ...
    
    
    glp_delete_prob(lp); // clean up
    delete []ia;
	delete []ja;
	delete []ar;
	
    return regret_ratio; 
}    
double worstDirection(int index, point_set_t *s, point_t* pt, float* &v)
{
	int K = index;
	int D = pt->dim;

    int* ia = new int[1+(K+5)*(D+5)];  //TODO: delete
	int* ja = new int[1+(K+5)*(D+5)];  //TODO: delete
    double* ar = new double[1+(K+5)*(D+5)];   //TODO: delete
    int i, j;   
    double epsilon=0.0000000000001; 


    glp_prob *lp;
    lp=glp_create_prob(); 
    glp_set_prob_name(lp, "max_regret_ratio");
    glp_set_obj_dir(lp, GLP_MAX);
    
 
    glp_add_rows(lp, K+2);  // add K+2 rows: q_1...q_k and r_1 and r_2
    // Add rows q_1 ... q_K 
    for (i=1; i<=K; i++) {
        char buf[10]; 
        sprintf(buf, "q%d", i);
        glp_set_row_name(lp, i, buf); 
        glp_set_row_bnds(lp, i, GLP_UP, 0.0, 0.0); // -infity < qi<=0
    }
       // Add rows r_1 and r_2
    glp_set_row_name(lp, K+1, "r1");  
    glp_set_row_bnds(lp, K+1, GLP_UP, 0.0, 1.0+epsilon); // r1 <= 1
    glp_set_row_name(lp, K+2, "r2");  
    glp_set_row_bnds(lp, K+2, GLP_UP, 0.0, -1.0+epsilon); // r2 <=-1
    
    
    glp_add_cols(lp, D+1);    // add D+1 columns: v[1] ... v[D] and x
    // Add col v[1] ... v[D]
    for (i=1; i<=D; i++) {
        char buf[10]; 
        sprintf(buf, "v%d", i);
        
        glp_set_col_name(lp, i, buf); 
        glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0); // 0 <= v[i] < infty
        glp_set_obj_coef(lp, i, 0.0);  // objective: maximize x ONLY
    }
    
    // Add col x
    glp_set_col_name(lp, D+1, "x");
    glp_set_col_bnds(lp, D+1, GLP_FR, 0.0, 0.0); // -infty <= x <= infty
    glp_set_obj_coef(lp, D+1, 1.0);  // objective: maximize x

    /*
            v[1]           v[2]   ...        v[D]    x
      q1 -(p-p_i)[0] -(p-p_i)[1] ... -(p-p_i)[D-1]   1
      q2  ......  same as q1  .................
      ...
      qK  ......  same as q1 .................
      r1   p[0]           p[1]    ...        p[D-1]  0
      r2  -p[0]         -p[1]    ...       -p[D-1]   0
    */
    
    int counter=1; 
    // set value on row q1 ... qk
    for (i=1; i<=K; i++) {
        
        #ifdef DEBUG_LP
        fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, i, D+1, 1.0); //DEBUG
        #endif
        
        ia[counter]=i; ja[counter]=D+1; ar[counter++]=1.0; // a["qi", "x"] =1
        for (j=1; j<=D; j++) {
            #ifdef DEBUG_LP
            fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, i, j, -(pt.a[j-1]-s[i-1].a[j-1])); //DEBUG
            #endif
            ia[counter]=i; ja[counter]=j; 
			ar[counter++] = -(pt->coord[j-1]-s->points[i-1]->coord[j-1]); // a["qi", "v[j]"] = -(pt-s[i-1])[j-1]
       }    
    }

    // set value on row r1 and r2
    #ifdef DEBUG_LP
    fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+1, D+1, 0); //DEBUG
    fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+2, D+1, 0); //DEBUG
    #endif
    ia[counter]=K+1; ja[counter]=D+1; ar[counter++]=0.0; // a["r1", "x"]=0 
    ia[counter]=K+2; ja[counter]=D+1; ar[counter++]=0.0; // a["r2", "x"]=0

    for (i=1; i<=D; i++){
        #ifdef DEBUG_LP
        fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+1, i, pt.a[i-1]); //DEBUG
        fprintf(stderr, "%d: a[%d, %d]=%lf\n", counter, K+2, i, -pt.a[i-1]); //DEBUG
        #endif
		ia[counter]=K+1; ja[counter]=i; ar[counter++]=pt->coord[i-1];   // e.g. a["r1", "v[1]"]=pt[0];   
		ia[counter]=K+2; ja[counter]=i; ar[counter++]=-pt->coord[i-1];   // e.g. a["r2", "v[1]"]=-pt[0];   
    }
    
    

    // loading data  
    glp_load_matrix(lp, counter-1, ia, ja, ar);    
    
    
    // Use this to print out the LP if you're debugging
    #ifdef DEBUG_LP
    glp_write_lp(lp, NULL, "testlp.lp");  // DEBUG
    #endif
    
    
    // running simplex
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev=GLP_MSG_OFF; // turn off all message by glp_simplex 
    glp_simplex(lp, &parm);    
    
    // write solution 
    #ifdef DEBUG_LP
    glp_print_sol(lp, "testlp.sol");  // DEBUG
    #endif
    
    // get values
    double regret_ratio=glp_get_obj_val(lp); 
    for (i=0; i<D; i++)
        v[i]=glp_get_col_prim(lp, i+1); // v[0] is at column 1, v[1] at col 2 ...
    
    
    glp_delete_prob(lp); // clean up
    delete []ia;
	delete []ja;
	delete []ar;
	
    return regret_ratio; 
}  

// determinant code from http://www.c.happycodings.com/Beginners_Lab_Assignments/code62.html
double determinant(int n, double** a)
{
    int i, j, k;
    double mult;
    double det = 1.0;

    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            mult = a[j][i]/a[i][i];
            for(k = 0; k < n; k++)
            {
                if (i==j) break;
                a[j][k] = a[j][k] - a[i][k] * mult;
            }
        }
    }

    for(i = 0; i < n; i++)
    {
        det = det * a[i][i];
    }

    return det;
}

void test1(){
    // This test should outputs regret ratio = 0.44444
    // and v= 0.5555 0.5555
 /*   struct point pt; 
    pt.d=2; 
    pt.a = (double*)malloc( 3*sizeof(double) ); 
    pt.a[0]=0.9;
    pt.a[1]=0.9; */

    point_t* pt = alloc_point(2);
	pt->coord[0] = 0.9;
	pt->coord[1] = 0.9;

	/*
    struct point s[2]; 
    s[0].d=2; 
    s[0].a=(double*)malloc( 3*sizeof(double) );
    s[0].a[0]=1; 
    s[0].a[1]=0; 
	*/

	point_set_t* s = alloc_point_set(2);
	s->points[0] = alloc_point(2);
	s->points[0]->coord[0] = 1;
	s->points[0]->coord[1] = 0;

	/*
    s[1].d=2; 
    s[1].a=(double*)malloc( 3*sizeof(double) );
    s[1].a[0]=0; 
    s[1].a[1]=1; */
    
	s->points[1] = alloc_point(2);
	s->points[1]->coord[0] = 0;
	s->points[1]->coord[1] = 1
		;
	
    double *v = new double[2]; 
    double regret_ratio=0;
    
    
    regret_ratio=worstDirection(s, pt, v); 
    printf("regret ratio = %lf\n", regret_ratio);    
    printf("v=%lf %lf\n", v[0], v[1]);  
}
void test2(){
    // This test should outputs negative regret ratio 
    // and v= 1.25 1.25
 /*   struct point pt; 
    pt.d=2; 
    pt.a = (double*)malloc( 3*sizeof(double) ); 
    pt.a[0]=0.4;
    pt.a[1]=0.4; 
 */
	point_t* pt = alloc_point(2);
	pt->coord[0] = 0.4;
	pt->coord[1] = 0.4;

/*
    struct point s[2]; 
    s[0].d=2; 
    s[0].a=(double*)malloc( 3*sizeof(double) );
    s[0].a[0]=1; 
    s[0].a[1]=0; 

    s[1].d=2; 
    s[1].a=(double*)malloc( 3*sizeof(double) );
    s[1].a[0]=0; 
    s[1].a[1]=1; 
    
    double v[2]; 
    double regret_ratio=0;
    
    
    regret_ratio=worstDirection(2, 2, s, pt, v); 
    printf("regret ratio = %lf\n", regret_ratio);    
    printf("v=%lf %lf\n", v[0], v[1]);
*/
	point_set_t* s = alloc_point_set(2);
	s->points[0] = alloc_point(2);
	s->points[0]->coord[0] = 1;
	s->points[0]->coord[1] = 0;

	s->points[1] = alloc_point(2);
	s->points[1]->coord[0] = 0;
	s->points[1]->coord[1] = 1
		;
	
    double *v = new double[2]; 
    double regret_ratio=0;
    
    
    regret_ratio=worstDirection(s, pt, v); 
    printf("regret ratio = %lf\n", regret_ratio);    
    printf("v=%lf %lf\n", v[0], v[1]); 
}

/*
* Compute the MRR of a given set of points
*/
double evaluateLP(point_set_t *p, point_set_t* S, int VERBOSE)
{
	int D = p->points[0]->dim;
	int N = p->numberOfPoints;
	int K = S->numberOfPoints;

	int i, j;
	double maxRegret = 0.0, maxK, maxN;
	double* v = new double[D];

	for (i = 0; i < N; ++i)
	{
		// obtain the worst utility vector
		worstDirection(S, p->points[i], v);

		maxN = dot_prod(maxPoint(p, v), v);
		maxK = dot_prod(maxPoint(S, v), v);

		if (1.0 - maxK / maxN > maxRegret)
			maxRegret = 1.0 - maxK / maxN;
	}

	if (VERBOSE)
		printf("LP max regret ratio = %lf\n", maxRegret);

	delete[]v;
	return maxRegret;
}

/*
int main() {
    test1();
    test2();        
    return 0; 
}
*/
