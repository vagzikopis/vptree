#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#ifndef VPTREE_H
#define VPTREE_H
#define ROWS 100
#define D 200
//vptree struct contains information about each
//node of the vantage point tree
struct vptree
{
  //vantage point's coordinates
  double *vp;
  //vantage point's index
  int idx;
  //vantage point's median
  double md;
  //outer subtree
  struct vptree* outer;
  //inner subtree
  struct vptree* inner;
};
typedef struct vptree vptree ;
/*_______________________________________________________________*/

//thread_distance contains arguments used in *calculator function
typedef struct{
  double *X;
  double *distances;
  int start;
  int stop;
  int n;
  int d;
}thread_distance;
/*_______________________________________________________________*/

//thread_tree struct contains arguments used in thread_function
typedef struct{
  int n;
  int d;
  int *index;
  double *X;
  vptree *T;
}thread_tree;
/*_______________________________________________________________*/

//This function calculates the distances from the vantage point using openMP framework
double *openmp_distances(double *X,int n,int d);
/*_______________________________________________________________*/

//This function calculates the distances from the vantage point sequentially
double * sequential_distances(double *X, int n, int d);
/*_______________________________________________________________*/

//This function calculates the distances from the vantage point using cilk framework
double * cilk_distances(double *X, int n, int d);
/*_______________________________________________________________*/

//This function calculates the distances from the vantage point using pthreads
//The number of threads created is based on n, the amount of distances
//that have to be calculated
double * pthread_distances(double *X, int n, int d);
/*_______________________________________________________________*/

//This function calculates the distances for each thread created
//by pthread_distances function
void *calculator(void *arg);
/*_______________________________________________________________*/

//This function is simply calls work function
//It is called when recursion is executed with pthreads
void * thread_function(void * args);
/*_______________________________________________________________*/

//This function returns an array with:
//array[0] equal to the number of elements with distances lower or equal with the median
//array[1] equal to the number of elements with distances greater than the median
int nInner(double *X, int nX, double median, double distances[]);
/*_______________________________________________________________*/

//This function is implements recursion and creates the vantage point tree
vptree *work(double *X, int n, int d, int index[]);
/*_______________________________________________________________*/

//The next three functions execute quickselect
double quickselect(double* arr, int l, int r, int k);
int partition(double* arr, int l, int r);
void swap(double* a, double* b);
/*_______________________________________________________________*/


//! Build vantage-point tree given input dataset X
/*!
\param X Input data points, stored as [n-by-d] array
\param n Number of data points (rows of X)
\param d Number of dimensions (columns of X)
\return The vantage-point tree
*/
vptree * buildvp(double *X, int n, int d);
//! Return vantage-point subtree with points inside radius
/*!
\param node A vantage-point tree
\return The vantage-point subtree
*/
vptree * getInner(vptree * T);
//! Return vantage-point subtree with points outside radius
/*!
\param node A vantage-point tree
\return The vantage-point subtree
*/
vptree * getOuter(vptree * T);
//! Return median of distances to vantage point
/*!
\param node A vantage-point tree
\return The median distance
*/
double getMD(vptree * T);
//! Return the coordinates of the vantage point
/*!
\param node A vantage-point tree
\return The coordinates [d-dimensional vector]
*/
double * getVP(vptree * T);
//! Return the index of the vantage point
/*!
\param node A vantage-point tree
\return The index to the input vector of data points
*/
int getIDX(vptree * T);
#endif
