# define the C/C++ compiler to use,default here is clang
CC = gcc-5

all:
	$(CC) -o ../vptree_sequential vptree_sequential.c main.c -lm; cd ..
	cd ..;./vptree_sequential
	$(CC) -o ../vptree_pthreads vptree_pthreads.c main.c -lm -pthread; cd ..
	cd ..;./vptree_pthreads
	$(CC) -o ../vptree_cilk vptree_cilk.c main.c -fcilkplus -lm; cd ..
	cd ..;./vptree_cilk
	$(CC) -o ../vptree_openmp vptree_openmp.c main.c -lm -fopenmp; cd ..
	cd ..;./vptree_openmp

	
