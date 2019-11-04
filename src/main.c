#include "../inc/vptree.h"

int main()
{
    int n, d;
    double *X, time;
    struct timeval startwtime, endwtime;

    n = 10000;
    d = 20;

    X = malloc(n*d*sizeof(double));

    //! Initialize data
    for (int i=0;i<n*d;i++)
      *(X + i) = (double)rand()/(double)RAND_MAX;

    gettimeofday (&startwtime, NULL);
    vptree *T;
    T = buildvp(X, n, d);
    gettimeofday (&endwtime, NULL);
    time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
  		      + endwtime.tv_sec - startwtime.tv_sec);

    printf("\t\tExecution Time: %f sec\n",time);

    free(X);
    return 0;
}
