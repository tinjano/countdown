/*
 * expr.h
 *
 *  Created on: Dec 12, 2023
 *      Author: tinjano
 *      https://github.com/tinjano/countdown
 *      git@github.com:tinjano/countdown.git
 */

#ifndef EXPR_H_
#define EXPR_H_

#include<stdlib.h>
#include<stdio.h>
#include<stdbool.h>
#include<omp.h>

#include<string.h>

#include<limits.h>
#include<sys/resource.h>

#define PLUS 43
#define MINUS 45
#define TIMES 42
#define DIV 47

#define START_SIZE_ 1000
#define START_SIZE_MAT_ 5

#define THRESHOLD_ ULLONG_MAX
// #define THRESHOLD_ 100000

#define MAX_MEMORY_ 9*1024*1024  // in KiB

typedef unsigned char op_t;
typedef unsigned long long value_t;
typedef unsigned char count_t;

/* vector */
typedef struct Vector{
	size_t size;
	count_t* val;
} Vector;

Vector* vector_init(size_t size) {
	Vector* vector = (Vector*)malloc(sizeof(Vector));

	if (vector == NULL) {
		fprintf(stderr, "Error: Memory allocation error (vector_init)\n");
		exit(1);
	}

	vector->val = (count_t*)malloc(size * sizeof(count_t));

	if (vector->val == NULL) {
			fprintf(stderr, "Error: Memory allocation error (vector_init)\n");
			exit(1);
	}

	memset(vector->val, 0, size * sizeof(count_t));
	vector->size = size;

	return vector;
}

void vector_sum(const Vector* a, const Vector* b, Vector* output) {
	if (a->size == b->size && a->size == output->size) {
		for (size_t i=0; i < a->size; ++i)
			output->val[i] = a->val[i] + b->val[i];
	} else {
		fprintf(stderr, "Error: vector sizes do not match (vector_sum)\n");
		exit(1);
	}
}

static inline int vector_isgt_all(const Vector* a, const Vector* b) {
	bool flag_gt = true, flag_eq = true;
	if (a->size == b->size) {
		for (size_t i=0; i < a->size; ++i) {
			flag_gt &= (a->val[i] >= b->val[i]);
			flag_eq &= (a->val[i] == b->val[i]);
		}
		return (int)flag_eq * 2 + (int)(flag_gt & !flag_eq);
	} else {
		fprintf(stderr, "Error: vector sizes do not match (vector_isgt)\n");
		exit(1);
	}
}

bool vector_isgt_any(const Vector* a, const Vector* b) {
	if (a->size == b->size) {
		bool flag = false;
		for (size_t i=0; i < a->size; ++i)
			flag |= (a->val[i] > b->val[i]);
		return flag;
	} else {
		fprintf(stderr, "Error: vector sizes do not match (vector_isgt)\n");
		exit(1);
	}
}

bool vector_sumgt(const Vector* a, const Vector* b, const Vector* c) {
	if (a -> size == b->size && c != NULL) {
		bool flag = true;
		for (size_t i=0; i < a->size; ++i)
			flag &= (a->val[i] + b->val[i] < c->val[i]);
		return flag;
	}
	else {
		fprintf(stderr, "Error: vector sizes do not match (vector_sumgt)\n");
		exit(1);
	}
}

count_t vector_elsum(const Vector* a) {
	count_t sum = 0;
	for (size_t i=0; i < a->size; ++i)
		sum += a->val[i];
	return sum;
}

void vector_print(const Vector* vector) {
	printf("Vector at address %p... ", vector);
	for (size_t i=0; i<vector->size; ++i)
		printf("%hu ", vector->val[i]);
	printf("\n");
}

void vector_free(Vector* vector) {
	if (vector != NULL) {
		free(vector->val);
		free(vector);
	}
}
/* vector */


/* map  */
typedef struct Map{
	struct Map* left;
	struct Map* right;

	value_t val;

	Vector** mat;
	size_t size;
	size_t capacity;
} Map;

Map* map_init(value_t val) {
	Map* map = (Map*)malloc(sizeof(Map));
	if (map == NULL) {
		fprintf(stderr, "Error: Memory allocation error (map_init)\n");
		exit(1);
	}

	map->left = NULL;
	map->right = NULL;
	map->val = val;
	map->size = 0;
	map->capacity = START_SIZE_MAT_;


	map->mat = (Vector**)malloc(START_SIZE_MAT_ * sizeof(Vector*));
	if (map->mat == NULL) {
		fprintf(stderr, "Error: Memory allocation error (map_init)\n");
		exit(1);
	}

	return map;
}

static inline bool map_insert_into_mat(Map* map, Vector* vector) {
	bool flag_subopt = false, flag_seen = false;
	for (size_t i=0; i < map->size; ++i) {
		int fs = vector_isgt_all(vector, map->mat[i]);
		flag_seen |= (fs/2);
		flag_subopt |= (fs%2);

		if (flag_subopt)
			break;
	}

	if (!flag_subopt && !flag_seen) {
		if (map->size == map->capacity) {
			map->capacity *= 20;
			map->mat = (Vector**)realloc(map->mat, map->capacity * sizeof(Vector*));

			if (map->mat == NULL){
				fprintf(stderr, "Error: Memory allocation error (map_insert_into_mat)\n");
				exit(1);
			}
		}

		map->mat[map->size++] = vector;
	}
	return !flag_subopt;
}

bool map_insert(Map* map, value_t val, Vector* vector) {
	while (1) {
		if (map->val == val)
			return map_insert_into_mat(map, vector);
		else if (val < map->val) {
			if (map->left != NULL)
				map = map->left;
			else {
				map->left = map_init(val);
				return map_insert_into_mat(map->left, vector);
			}
		} else {
			if (map->right != NULL)
				map = map->right;
			else {
				map->right = map_init(val);
				return map_insert_into_mat(map->right, vector);
			}
		}
	}
}

void map_print(Map* map) {
	if (map == NULL)
		return;
	else {
		printf("Map at address %p...\n", map);
		printf("%llu has vectors:\n", map->val);
		for (int i=0; i < map->size; ++i)
			vector_print(map->mat[i]);
		printf("\n\n");
		map_print(map->left);
		map_print(map->right);
	}
}



void map_free(Map* map) {
	if (map != NULL) {
//		for (size_t i=0; i < map->size; ++i)
//			vector_free(map->mat[i]);
		free(map->mat);
		map_free(map->left);
		map_free(map->right);
		free(map);
	}
}
/* map */


/* expression */
typedef struct ExprTree{
	value_t val;
	Vector* vector;
	struct ExprTree* left;
	op_t o;
	struct ExprTree* right;
} ExprTree;

ExprTree* exprtree_init(value_t val, Vector* vector) {
	ExprTree* expr = (ExprTree*)malloc(sizeof(ExprTree));
	if (expr == NULL) {
		fprintf(stderr, "Error: Memory allocation error (ExprTree)\n");
		exit(1);
	}

	expr->val = val;
	expr->vector = vector;
	expr->left = NULL;
	expr->right = NULL;
	expr->o = 0;

	return expr;
}

void exprtree_printexpr(const ExprTree* expr) {
	if (expr->left == NULL)
		printf("%llu", expr->val);
	else {
		if (expr->o == TIMES || expr->o == DIV) {
			if (expr->left->o == PLUS || expr->left->o == MINUS) {
				printf("(");
				exprtree_printexpr(expr->left);
				printf(")");
			} else exprtree_printexpr(expr->left);
			printf("%c", expr->o);
			if (expr->right->o == PLUS || expr->right->o == MINUS) {
				printf("(");
				exprtree_printexpr(expr->right);
				printf(")");
			} else exprtree_printexpr(expr->right);
		} else {
			exprtree_printexpr(expr->left);
			printf("%c", expr->o);
			exprtree_printexpr(expr->right);
		}
	}
}

void exprtree_free(ExprTree* expr) {
	if (expr != NULL) {
		vector_free(expr->vector);
//		exprtree_free(expr->left);
//		exprtree_free(expr->right);
		free(expr);
	}
}
/* expression */


/* forest */
typedef struct ExprForest{
	ExprTree** expr;
	Map* map;
	Vector* maxvec;
	size_t* jumps;
	size_t top;
	size_t capacity;
} ExprForest;

ExprForest* exprforest_init() {
	ExprForest* forest = (ExprForest*)malloc(sizeof(ExprForest));
	if (forest == NULL) {
		fprintf(stderr, "Error: Memory allocation error (exprforest)\n");
		exit(1);
	}

	forest->expr = (ExprTree**)malloc(START_SIZE_ * sizeof(ExprTree*));
	if (forest->expr == NULL) {
		fprintf(stderr, "Error: Memory allocation error(exprforest expr)\n");
		exit(1);
	}

	forest->map = NULL;
	forest->maxvec = NULL;
	forest->jumps = NULL;
	forest->top = 0;
	forest->capacity = START_SIZE_;

	return forest;
}

void exprforest_free(ExprForest*);
void exprforest_print(ExprForest*);

void exprforest_insert(ExprForest* forest, ExprTree* expr) {
	if (forest->top == forest->capacity) {
		forest->capacity *= 50;
		forest->expr = (ExprTree**)realloc(forest->expr, forest->capacity * sizeof(ExprTree*));

		if (forest->expr == NULL) {
			fprintf(stderr, "Error: Memory Reallocation error (exprforest_insert)\n");
			exit(1);
		}
	}

	forest->expr[forest->top++] = expr;
}

void exprforest_input(ExprForest* forest, int N, int args[]) {
	forest->maxvec = vector_init(N);

	size_t u = 0;
	count_t temp[N];

	for (size_t i=0; i<N; ++i) {
		bool unique_flag = 1;
		int place = -1;
		for (size_t j=0; j<u; ++j)
			if (args[i] == temp[j]) {
				unique_flag = 0;
				place = j;
			}

		if (unique_flag) {
			Vector *vector = vector_init(N);
			vector->val[u] = 1;

			ExprTree* expr = exprtree_init(args[i], vector);
			exprforest_insert(forest, expr);

			forest->maxvec->val[u]++;
			temp[u++] = args[i];

		} else
			forest->maxvec->val[place]++;
	}

	forest->map = map_init(0);

	forest->maxvec->size = u;
	for (size_t i=0; i < forest->top; ++i) {
		forest->expr[i]->vector->size = u;
		map_insert(forest->map, forest->expr[i]->val, forest->expr[i]->vector);
	}

	forest->jumps = (size_t*)malloc((N+1) * sizeof(size_t));
	if (forest->jumps == NULL) {
		fprintf(stderr, "Error: Memory allocation error (jumps)\n");
		exit(1);
	}
	forest->jumps[1] = 0;
	forest->jumps[2] = u;
}

void exprforest_grow(ExprForest* forest, count_t target_sum) {
	size_t top = forest->top;
	for (size_t i=0; i<top; ++i) {

		ExprTree* first = forest->expr[i];
		count_t i_sum = vector_elsum(first->vector);
		size_t L = forest->jumps[target_sum - i_sum], U = forest->jumps[target_sum - i_sum + 1];

		for (size_t j=L; j<U; ++j) {
			ExprTree* second = forest->expr[j];

			if (first->val < second->val || (first->val == second->val && first < second))
				continue;

			Vector* combvec = vector_init(first->vector->size);
			vector_sum(first->vector, second->vector, combvec);

			if (vector_isgt_any(combvec, forest->maxvec))
				continue;


			bool decreasing = first->right == NULL
			                || first->right->val > second -> val
			                || (first->right->val == second->val && first->right > second);

			op_t o = PLUS;
			if (first->val + second->val < THRESHOLD_ && second->o != PLUS && second->o != MINUS
			   && ((first->o == TIMES || first->o == DIV || first->o == 0)
			       || (first->o == PLUS && decreasing))) {

				if (map_insert(forest->map, first->val + second->val, combvec)) {
					ExprTree* new_expr = exprtree_init(first->val + second->val, combvec);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					exprforest_insert(forest, new_expr);
				}
			}


			o = MINUS;
			if (first->val > second->val && second->o != PLUS && second->o != MINUS
			   && ((first->o == PLUS || first->o == TIMES || first->o == DIV || first->o == 0)
			       || (first->o == MINUS && decreasing))) {

				if (map_insert(forest->map, first->val - second->val, combvec)) {
					ExprTree* new_expr = exprtree_init(first->val - second->val, combvec);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					exprforest_insert(forest, new_expr);
				}
			}


			o = TIMES;
			if (second->o != TIMES && second->o != DIV
			   && ((first->o == PLUS || first->o == MINUS || first->o == 0)
			       || (first->o == TIMES && decreasing))) {

				if (first->val * second->val < THRESHOLD_ && map_insert(forest->map, first->val * second->val, combvec)) {
					ExprTree* new_expr = exprtree_init(first->val * second->val, combvec);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					exprforest_insert(forest, new_expr);
				}
			}


			o = DIV;
			if (first->val % second->val == 0 && second->o != TIMES && second->o != DIV
			   && ((first->o == PLUS || first->o == MINUS || first->o == TIMES || first->o == 0)
			       || (first->o == DIV && decreasing))) {

				if (map_insert(forest->map, first->val / second->val, combvec)) {
					ExprTree* new_expr = exprtree_init(first->val / second->val, combvec);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					exprforest_insert(forest, new_expr);
				}
			}
		}

		if (i%1000 == 0) {
			struct rusage usage;
			getrusage(RUSAGE_SELF, &usage);
			if (usage.ru_maxrss > MAX_MEMORY_) {
				fprintf(stderr, "Memory limit exceeded (%d MiB). The program will terminate.\n", MAX_MEMORY_/1024);
				exprforest_free(forest);
				exit(1);
			}
		}
	}
	forest->jumps[target_sum + 1] = forest->top;
}

void exprforest_query(const ExprForest* forest, value_t query) {
	int count = 0;
	printf("------\nQuery: %llu\n", query);
	for (size_t i=0; i < forest->top; ++i)
		if (forest->expr[i]->val == query) {
			exprtree_printexpr(forest->expr[i]);
			++count;
			printf("\n");
		}
	printf("%d solutions found.\n------\n", count);
}

void exprforest_print(ExprForest* forest) {
	printf("Exprforest at address %p... top %zu capacity %zu\n", forest, forest->top, forest->capacity);

	for (size_t i=0; i<forest->top; ++i) {
		printf("\tExpression with value %llu\n", forest->expr[i]->val);
		printf("\tHas ");
		vector_print(forest->expr[i]->vector);
		printf("\t");
		exprtree_printexpr(forest->expr[i]);
		printf("\n\n");
	}
}

void exprforest_free(ExprForest* forest) {
	for (size_t i=0; i < forest->top; ++i)
		forest->expr[i]->vector->val[0] = i;
	for (size_t i=0; i < forest->top; ++i)
		if (i == forest->expr[i]->vector->val[0])
			vector_free(forest->expr[i]->vector);

	free(forest->expr);
	map_free(forest->map);
	vector_free(forest->maxvec);
	free(forest);
}
/* forest */


#endif /* EXPR_H_ */
