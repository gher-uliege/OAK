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

Test Simple wrapper for CHOLMOD

*/
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#define INDEX size_t
#define N 4


int solve(INDEX n,INDEX nz, int* Si, int* Sj,double* Ss, double* bb, double* xx);

int main(void)
{
    int i;
    size_t nnz, nz;

    int* Si;
    int* Sj;
    double* Ss;
    double* bb;
    double* xx;

    nz = N;

    Si = (int*)malloc(sizeof(int)*nz);
    Sj = (int*)malloc(sizeof(int)*nz);
    Ss = (double*)malloc(sizeof(double)*nz);
    bb = (double*)malloc(sizeof(double)*N);
    xx = (double*)malloc(sizeof(double)*N);

    for (i = 0; i < nz; i++) {
      Si[i] = i;
      Sj[i] = i;
      Ss[i] = i+1;
    }

    for (i = 0; i < N; i++) {
      bb[i] = i;
    }
    
    solve(N,nz,Si,Sj,Ss,bb,xx);

    for (i = 0; i < N; i++) {
      printf ("x  %g\n",xx[i]);
    }

    free(Si);
    free(Sj);
    free(Ss);
    free(bb);
    free(xx);


    return (0);
}
