/*
!
!  OAK, Ocean Assimilation Kit
!  Copyright(c) 2015 Alexander Barth and Luc Vandenblucke
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!

Simple wrapper for CHOLMOD

 */

#include "cholmod.h"
#define INDEX size_t

/* we are working with pointer to pointer so that the Fortran code does not need 
   to know about the cholmod structures and that the C-code can modify the pointers*/

int cholmod_wrapper_start(cholmod_common **c) {  
  *c = malloc(sizeof(cholmod_common));
  cholmod_start(*c);
  /*printf("here %d\n",c);
    printf("here %d\n",*c);*/

  return 0;
}

int cholmod_wrapper_finish(cholmod_common **c) {  
  /*printf("almost done %d\n",*c);
    printf("done\n");*/
  cholmod_finish(*c);
  free(*c);

  return 0;
}


int cholmod_wrapper_matrix(cholmod_common **c,INDEX n,INDEX nz, int* Si, int* Sj,double* Ss, cholmod_sparse **A) {
    cholmod_triplet *At;    
    cholmod_dense *x, *b;
    cholmod_factor *L;
    int i;

    // nrow,ncol,nzmax,storage type (1=upper)
    At = cholmod_allocate_triplet(n, n, nz, 1, CHOLMOD_REAL, *c);

    At->nnz = nz;
    for (i = 0; i < nz; i++) {
      ((int*)At->i)[i] = Si[i];
      ((int*)At->j)[i] = Sj[i];
      ((double*)At->x)[i] = Ss[i];
    }

    *A = cholmod_triplet_to_sparse(At,At->nzmax,*c);
    cholmod_free_triplet(&At, *c);

    return 0;
}

int cholmod_wrapper_factorize(cholmod_common **c,cholmod_sparse **A, cholmod_factor **L) {

    *L = cholmod_analyze(*A, *c);		    /* analyze */
    cholmod_factorize(*A, *L, *c);		    /* factorize */

    return (*c)->status;
}


int cholmod_wrapper_solve(cholmod_common **c,cholmod_sparse **A, cholmod_factor **L, double* bb, double* xx) {
    cholmod_dense *x, *b;
    int i;

    b = cholmod_zeros((*A)->nrow, 1, (*A)->xtype, *c);
    for (i = 0; i < (*A)->nrow; i++) {
      ((double*)b->x)[i] = bb[i];
    }

    x = cholmod_solve(CHOLMOD_A, *L, b, *c);	    /* solve Ax=b */

    for (i = 0; i < (*A)->nrow; i++) {
      xx[i] = ((double*)x->x)[i];
    }

    cholmod_free_dense(&x, *c);
    cholmod_free_dense(&b, *c);

    return 0;
}

int cholmod_wrapper_free(cholmod_common **c,cholmod_sparse **A, cholmod_factor **L) {
    cholmod_free_factor(L, *c);
    cholmod_free_sparse(A, *c);

    return 0;
}


/*
solves A x = b for x where A is a symetric positive defined matrix

Input
n: size of A
Si: i-index of non-zeros elements in A (upper part)
Sj: j-index of non-zeros elements in A (upper part)
Ss: values of non-zeros elements in A (upper part)
bb: vector b

Output:
xx: vector x
*/


int solve(INDEX n,INDEX nz, int* Si, int* Sj,double* Ss, double* bb, double* xx) {
    cholmod_sparse *A;
    cholmod_triplet *At;    
    cholmod_dense *x, *b;
    cholmod_factor *L;
#ifdef DEBUG
    /* variables to check residual */
    cholmod_dense *r;
    double one [2] = {1,0}, m1 [2] = {-1,0};
#endif
    cholmod_common c;
    int i;

    /* init CHOLMOD */
    cholmod_start(&c);

    // nrow,ncol,nzmax,storage type (1=upper)
    At = cholmod_allocate_triplet(n, n, nz, 1, CHOLMOD_REAL, &c);

    At->nnz = nz;
    for (i = 0; i < nz; i++) {
      ((int*)At->i)[i] = Si[i];
      ((int*)At->j)[i] = Sj[i];
      ((double*)At->x)[i] = Ss[i];
    }

    A = cholmod_triplet_to_sparse(At,At->nzmax,&c);
    cholmod_free_triplet(&At, &c);

    /* A must be symmetric */
    if (A == NULL || A->stype == 0)
    {
        cholmod_free_sparse(&A, &c);
	cholmod_finish(&c);
	return -1;
    }

    b = cholmod_zeros(A->nrow, 1, A->xtype, &c);
    for (i = 0; i < n; i++) {
      ((double*)b->x)[i] = bb[i];
    }

    L = cholmod_analyze(A, &c);		    /* analyze */
    cholmod_factorize(A, L, &c);		    /* factorize */
    x = cholmod_solve(CHOLMOD_A, L, b, &c);	    /* solve Ax=b */

    for (i = 0; i < n; i++) {
      xx[i] = ((double*)x->x)[i];
    }

#ifdef DEBUG
    r = cholmod_copy_dense(b, &c);		    /* r = b */
    cholmod_sdmult(A, 0, m1, one, x, r, &c);	    /* r = r-Ax */
    printf("norm(b-Ax) %8.1e\n",
	    cholmod_norm_dense(r, 0, &c));	    /* print norm(r) */
    cholmod_free_dense(&r, &c);
#endif

    /* free matrices and vectors */
    cholmod_free_factor(&L, &c);
    cholmod_free_sparse(&A, &c);
    cholmod_free_dense(&x, &c);
    cholmod_free_dense(&b, &c);
    cholmod_finish(&c);
    return 0;
}

