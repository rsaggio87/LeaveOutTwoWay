/* ========================================================================== */
/* === CMG/Include/cmg.h          =========================================== */
/* ========================================================================== */

/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributd under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */



/* 	
 *  CMG base data types
 */

typedef unsigned int mIndex;
typedef unsigned int mSize;
typedef unsigned int boolean;



typedef struct smatrix {
    float    *a;
    mIndex   *ia; 
    mIndex   *ja;
    mSize    n;
    boolean  issym;       /* true if matrix is symmetric */
} s_matrix ;


typedef struct dmatrix {
    double    *a;
    mIndex   *ia; 
    mIndex   *ja;
    mSize    n;
    boolean  issym;       /* true if matrix is symmetric */
} d_matrix ;

#ifdef SINGLE_PR
typedef float precision;
typedef s_matrix matrix;
#else
typedef double precision;
typedef d_matrix matrix;
#endif


typedef struct cholesky {
    d_matrix ld;
    d_matrix ldT;
    mSize    *p;
    mSize    *invp;
} ldl_p ;

typedef struct shlevel {
    mIndex    *cI;            /* cluster index */ 
    mSize     nc;             /* number of clusters */
    matrix    A;              /* single precision symmetric matrix */
    precision *invD;          /* inverse of diagonal  */
    boolean   dc;             /* true if diagonal correction occurs */
    int       repeat;         /* number of recursive calls to next level */

    /* add update info */
    boolean       islast;     /* true for last level */
    boolean       laplacian;  /* true if system is laplacian (1st level only) */
    boolean       iterative;  /* true if an iterative method is used */
    ldl_p         chol;       /* cholesky factorization   - last level only */
    
    precision     *lws1,*lws2,*sws1,*sws2, *sws3;   /* working space for the solver */
} s_hlevel ; 


/* entry-wise vector operations */
void vmv   (precision *, precision *, precision *, mSize);
void vpv   (precision *, precision *, precision *, mSize);
void vvmul (precision *, precision *, precision *, mSize);
void vpvmv (precision *, precision *, precision *, precision *, mSize);

/* matrix-vector multiplication prototype */
void spmv (mSize, precision *, mIndex *, mIndex *, precision *, precision *);
void sspmv (mSize, precision *, mIndex *, mIndex *, precision *, precision *);


/* restriction operators */
void rmvecmul (mIndex *, precision *, mSize ,  precision *, mSize , boolean);
void trmvecmul (mIndex *, precision *, mSize , precision *, mSize , boolean);

/* factorized solve */
void ldl_solve (ldl_p *, precision *, precision *);

/* main solver */
void preconditioner ( s_hlevel *, precision *, int , int, precision *);


