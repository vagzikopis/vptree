// Sequential creation of a Vantage-Point Tree
// Author : Zikopis Evangelos  <vagzikopis@gmail.com>
// Date : 3/11/2019
// Function documentation at vptree.h

#include "../inc/vptree.h"
#include <pthread.h>

//Variable counting active threads
int active_threads = 0;
//Maximum thread threshold
int max_threads = 10;


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
  //Find and store the distances from Vantage Point using pthreads
  //If n>1000 calculate distances sequentially
  double *distances;
  if (n > 10000 )
  {
    distances = pthread_distances(X,n,d);
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

  //Control active threads with this if clause
  if (n > 100 && active_threads < max_threads)
  {
    pthread_t pid;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    //Create a struct in order to pass arguments to the thread function
    thread_tree *arg = calloc(1,sizeof(thread_tree));
    arg->n = innerLength;
    arg->d = d;
    arg->X = Inner;
    arg->index = indexIn;

    pthread_t pid_2;
    pthread_attr_t attr_2;
    pthread_attr_init(&attr_2);
    //Create a struct in order to pass arguments to the thread function
    thread_tree *arg_2 = calloc(1,sizeof(thread_tree));
    arg_2->n = outerLength;
    arg_2->d = d;
    arg_2->X = Outer;
    arg_2->index = indexOut;
    //Create the first thread in order to calculate the inner subtree
    int rc = pthread_create(&pid,&attr,thread_function,arg);
    if(rc)
    {
      printf("Error pthread_create: %d\nhere\n", rc);
      exit(-1);
    }
    //Create the second thread in order to calculate the outer subtree
    int rc_2 = pthread_create(&pid_2,&attr_2,thread_function,arg_2);
    if(rc_2)
    {
      printf("Error pthread_create: %d\nhere\n", rc);
      exit(-1);
    }
    //Wait threads
    pthread_join(pid,NULL);
    pthread_join(pid_2,NULL);
    T->inner = arg->T;
    T->outer = arg_2->T;
    //Free allocated memory
    free(arg);
    free(arg_2);

  }else
  {
    //If active_threads >= max_threads then work recursively
    T->inner = work(Inner,innerLength,d,indexIn);
    T->outer = work(Outer,outerLength,d,indexOut);
  }
  //Free allocated memory
  free(Inner);
  free(Outer);
  free(indexIn);
  free(indexOut);
  free(distances);
  return T;
}
void * thread_function(void * args)
{
  thread_tree *arg = (thread_tree *) args;
  //Increase the number of active threads
  active_threads++;
  arg->T = work(arg->X,arg->n,arg->d,arg->index);
  //Decrease the number of active threads
  active_threads--;
  pthread_exit(NULL);
}

void *calculator(void *arg)
{
  thread_distance *arg_struct = (thread_distance *) arg;
  int n = arg_struct->n;
  double *X = arg_struct->X;
  double *distances = arg_struct->distances;
  int start = arg_struct->start;
  int stop = arg_struct->stop;
  int d = arg_struct->d;
  for (int i=start; i<stop; i++)
  {
    for (int j=0; j<d; j++)
    {
      *(distances+i) += pow(*(X+n*j+n-1) - *(X+n*j+i), 2);
    }
    *(distances+i) = sqrt(*(distances+i));
  }
  pthread_exit(NULL);
}

double * pthread_distances(double *X, int n, int d)
{
  double *distances = calloc(n-1,sizeof(double));
  int thread_num;
  //Calculate thread_num according to the number of points
  //The ideal thread_num depends on machine's hardware
  if (n-1 > 100000)
  {
    thread_num = 6;
  }else if (n-1 > 50000)
  {
    thread_num = 4;
  }else if (n-1 > 10000)
  {
    thread_num = 2;
  }else
  {
    thread_num = 1;
  }
  pthread_t pids[thread_num];
  thread_distance *args=malloc(thread_num * sizeof(thread_distance));
  //Create structs to pass arguments to *calculator function
  for (int i=0; i<thread_num; i++)
  {
    args[i].X = X;
    args[i].n = n;
    args[i].d = d;
    //Start and Stop inform the thread about the memory cells
    //that the result should be stored
    args[i].start = (int)(((n-1)/thread_num)*i);
    args[i].stop = args[i].start + (int)((n-1)/thread_num);
    args[i].distances = distances;
    if ((n-1)%thread_num != 0 && i==thread_num-1)
    {
      args[i].stop = n-1;
    }
  }
  int rc;
  //Create threads
  for (int i=0; i<thread_num; i++)
  {
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    rc = pthread_create(&pids[i], &attr, calculator, &args[i]);
    if (rc)
    {
      printf("Pthread error code: %d \n ", rc);
      exit(-1);
    }
  }
  //Wait for threads to join
  for (int i=0; i<thread_num; i++)
  {
    pthread_join(pids[i], NULL);
  }
  free(args);
  return distances;
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
