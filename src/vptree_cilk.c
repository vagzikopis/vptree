// Sequential creation of a Vantage-Point Tree
// Author : Zikopis Evangelos  <vagzikopis@gmail.com>
// Date : 3/11/2019
// Function documentation at vptree.h

#include "../inc/vptree.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>


vptree * buildvp(double *X, int n, int d)
{
  //index stores the indexes of all points
  int *index;
  index = malloc(n*sizeof(int));
  for (int i=0; i<n; i++)
  {
    index[i] = i;
  }
  vptree *T = work(X,n,d,index);
  free(index);
  return T;
}

vptree * work(double *X, int n, int d, int index[])
{
  vptree *T = NULL;
  //If the number of points is zero return an empty vptree
  if ( n == 0)
  {
    return T;
  }
  //Allocate memory
  T = calloc(1,sizeof(vptree));
  T->vp = malloc(d*sizeof(double));
  //Store Vantage Point's coordinates
  for (int i=0; i<d; i++)
  {
    T->vp[i] = *(X + n-1 + i*n);
  }
  //Store Vantage Point's index
  T->idx = index[n-1];
  if ( n == 1)
  {
    return T;
  }
  //Find and store the distances from Vantage Point
  //If n > 1000 calculate distances using Cilk
  //Otherwise sequentially
  double *distances;
  if(n>1000)
  {
    distances = cilk_distances(X,n,d);
  }else
  {
    distances = sequential_distances(X,n,d);
  }
  double *temp = malloc((n-1)*sizeof(double));
  for(int i=0; i<n-1; i++){
    temp[i] = distances[i];
  }
  //Select the median distance with quickselect
  if ( (n-1)%2 == 0 )
  {
    T->md = quickselect(temp,0,n-2,(n-1)/2+1);
    T->md += quickselect(temp,0,n-2,(n-1)/2);
    T->md /= 2.0;
  }else
  {
    T->md = quickselect(temp,0,n-2,n/2);
  }
  free(temp);
  //Find the number of inner and outer points
  int innerLength = nInner(X,n,T->md,distances);
  int outerLength = n-1 -innerLength;
  //Allocate memory for the two new sets
  double *Outer = malloc(outerLength*d*sizeof(double));
  double *Inner = malloc(innerLength*d*sizeof(double));
  int *indexIn = malloc(innerLength*sizeof(int));
  int *indexOut = malloc(outerLength*sizeof(int));
  int innerCounter = 0;
  int outerCounter = 0;
  //Fullfill the arrays with the appropriate points
  for(int i=0; i<n-1; i++)
  {
    if(distances[i] <= T->md)
    {
      indexIn[innerCounter] = index[i];
      for(int j=0; j<d; j++)
      {
        Inner[innerCounter + innerLength*j] = X[n*j+i];
      }
      innerCounter++;
    }else
    {
      indexOut[outerCounter] = index[i];
      for(int j=0; j<d; j++)
      {
        Outer[outerCounter + outerLength*j] = X[n*j+i];
      }
      outerCounter++;
    }
  }
  //Use cilk_spawn
  T->inner = cilk_spawn work(Inner,innerLength,d,indexIn);
  //Calculate the outer subtree recursively
  T->outer = work(Outer,outerLength,d,indexOut);
  cilk_sync;
  //Free allocated memory
  free(Inner);
  free(Outer);
  free(indexIn);
  free(indexOut);
  free(distances);
  return T;
}

int nInner(double *X, int nX, double median, double distances[])
{
  int result, nInner;
  nInner = 0;

  for (int i=0; i<nX-1; i++)
  {
    if (distances[i] <= median)
    {
      nInner++;
    }
  }
  result = nInner;
  return result;
}

double * cilk_distances(double *X, int n, int d)
{
  double *distances = calloc(n-1,sizeof(double));
  //Calculate (xvp - xi)^2, (yvp - yi)^2, ... and then square the summary
  //Parallelise the for loop with cilk_for
  cilk_for (int i=0; i<n-1; i++)
  {
    for (int j=0; j<d; j++)
    {
      distances[i] += pow(*(X+n*j+n-1) - *(X+n*j+i), 2);
    }
    distances[i] = sqrt(distances[i]);
  }
  return distances;
}

double * sequential_distances(double *X, int n, int d)
{
  double *distances = calloc(n-1,sizeof(double));
  //Calculate (xvp - xi)^2, (yvp - yi)^2, ... and then square the summary
  for (int i=0; i<n-1; i++)
  {
    for (int j=0; j<d; j++)
    {
      distances[i] += pow(*(X+n*j+n-1) - *(X+n*j+i), 2);
    }
    distances[i] = sqrt(distances[i]);
  }
  return distances;
}

vptree * getInner(vptree *T)
{
  return T->inner;
}
vptree * getOuter(vptree *T)
{
  return T->outer;
}
double getMD(vptree *T)
{
  return T->md;
}
double * getVP(vptree *T)
{
  return T->vp;
}
int getIDX(vptree *T)
{
  return T->idx;
}

void swap(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

int partition(double* arr, int l, int r)
{
    double x = *(arr + r);
    int i = l;
    for (int j = l; j <= r - 1; j++)
    {
        if (*(arr + j) <= x)
        {
            swap(&*(arr + i), &*(arr + j));
            i++;
        }
    }
    swap(&*(arr + i), &*(arr + r));
    return i;
}

double quickselect(double* arr, int l, int r, int k)
{
    if (k > 0 && k <= r - l + 1)
    {
        int index = partition(arr, l, r);

        if (index - l == k - 1)
            return *(arr + index);

        if (index - l > k - 1)
            return quickselect(arr, l, index - 1, k);

        return quickselect(arr, index + 1, r, k - index + l - 1);
    }

    return 0;
}
