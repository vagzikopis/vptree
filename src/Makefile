# define the C/C++ compiler to use,default here is clang
CC = gcc-5

all:
	cd src; $(CC) -o ../vptree_sequential vptree_sequential.c main.c -lm; cd ..
	./vptree_sequential
	cd src; $(CC) -o ../vptree_pthreads vptree_pthreads.c main.c -lm -pthread; cd ..
	./vptree_pthreads
	cd src; $(CC) -o ../vptree_cilk vptree_cilk.c main.c -fcilkplus -lm; cd ..
	./vptree_cilk
	cd src; $(CC) -o ../vptree_openmp vptree_openmp.c main.c -lm -fopenmp; cd ..
	./vptree_openmp

	
