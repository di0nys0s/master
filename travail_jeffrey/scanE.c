#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
int main (int argc, char *argv[]) {
  int th_id, nthreads,n,chunk;
  int npos=2048;
  double w=0.008596858979444764;                          
  double xmin=0.02;
  double xmax=30.;
  double x[npos];
  double delr;
  chunk=8;
  delr=(xmax-xmin)/(npos-1);
  #pragma omp parallel shared(x,nthreads) private(n,th_id)
{
 th_id = omp_get_thread_num();
  if (th_id == 0)
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d starting...\n",th_id);




        for ( n = 1; n <= npos; n++ ) {
          x[n]=xmin+(n-1)*delr;
          printf ("x(%d) vaut: %f \n", n,x[n]);
        }


    }  /* end of sections */

  return 0;

}
