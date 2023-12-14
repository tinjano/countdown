/*
 * countdown.c
 *
 *  Created on: Dec 8, 2023
 *      Author: tinjano
 *      https://github.com/tinjano/countdown
 *      git@github.com:tinjano/countdown.git
 */

#include<stdio.h>
#include<time.h>
#include "expr.h"

int main(int argc, char* argv[]){
	size_t N, Q;

	if (argc < 5) {
		fprintf(stderr, "Usage %s num_numbers num_queries <numbers...> <queries...>\n", argv[0]);
		return 1;
	}

	N = atoi(argv[1]);
	Q = atoi(argv[2]);

	if (N < 1 || Q < 0) {
		fprintf(stderr, "Invalid inputs for N and Q\n");
		return 1;
	}

	if (argc != 3 + N + Q) {
	    fprintf(stderr, "Error: Incorrect number of integers provided\n");
	    return 1;
	}

	int numbers[N], q;

	ExprForest* forest = exprforest_init();

	for (size_t i=0; i<N; ++i) {
		numbers[i] = atoi(argv[i+3]);
		if (numbers[i] < 1) {
			fprintf(stderr, "Error: Numbers must be positive integers\n");
			return 1;
		}
	}

	exprforest_input(forest, N, numbers);

	clock_t start = clock();
	for (int i=0; i<N-1; ++i)
		exprforest_grow(forest, i+2);

	for (size_t i=0; i<Q; ++i) {
		q = atoi(argv[3+N+i]);
		exprforest_query(forest, q);
	}

	printf("%zu solutions found in total.\n", forest->top);
	clock_t end = clock();
	double time = (double)(end-start)/CLOCKS_PER_SEC;
	printf("Time:%lf seconds", time);
	if (time < 1.0)
		printf(" (%lf miliseconds)\n", time*1000.0);
	else
		printf("\n");

	exprforest_free(forest);
	printf("Memory freed successfully.\n");
	return 0;
}

