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
#include<stdint.h>
#include<stdbool.h>
#include<omp.h>

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
typedef unsigned long long mask_t;

/* map  */
typedef struct Map{
	struct Map* left;
	struct Map* right;

	value_t val;

	mask_t* masks;
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


	map->masks = (mask_t*)malloc(START_SIZE_MAT_ * sizeof(mask_t));
	if (map->masks == NULL) {
		fprintf(stderr, "Error: Memory allocation error (map_init)\n");
		exit(1);
	}

	return map;
}

static inline bool map_insert_into_mat(Map* map, mask_t mask) {
	bool flag_subopt = false, flag_seen = false, diff;
	mask_t* masks = map->masks;
	for (size_t i=0; i < map->size; ++i) {
		diff = mask != masks[i];

		flag_subopt |= (diff && mask % masks[i] == 0);
		flag_seen |= (!diff);

		if (flag_subopt)
			break;
	}

	if (!flag_subopt && !flag_seen) {
		if (map->size == map->capacity) {
			map->capacity *= 20;
			mask_t* temp = (mask_t*)realloc(map->masks, map->capacity * sizeof(mask_t));

			if (temp == NULL){
				free(map->masks);
				fprintf(stderr, "Error: Memory allocation error (map_insert_into_mat)\n");
				exit(1);
			}

			map->masks = temp;
		}

		map->masks[map->size++] = mask;
	}
	return !flag_subopt;
}

bool map_insert(Map* map, value_t val, mask_t mask) {
	while (1) {
		if (map->val == val)
			return map_insert_into_mat(map, mask);
		else if (val < map->val) {
			if (map->left != NULL)
				map = map->left;
			else {
				map->left = map_init(val);
				return map_insert_into_mat(map->left, mask);
			}
		} else {
			if (map->right != NULL)
				map = map->right;
			else {
				map->right = map_init(val);
				return map_insert_into_mat(map->right, mask);
			}
		}
	}
}

int map_heights(Map* map, int depth) {
	if (map == NULL)
		return -1;
	size_t lh = 1 + map_heights(map->left, depth+1);
	size_t rh = 1 + map_heights(map->right, depth+1);

	for (int i=0; i<depth; ++i)
		printf("\t");
	printf("Address %p value %llu -- left: %lu, right: %lu\n", map, map->val, lh, rh);
	return lh + rh;
}

void map_print(Map* map) {
	if (map == NULL)
		return;
	else {
		printf("Map at address %p...\n", map);
		printf("%llu has masks:\n", map->val);
		for (int i=0; i < map->size; ++i)
			printf("%llu ", map->masks[i]);
		printf("\n\n");
		map_print(map->left);
		map_print(map->right);
	}
}



void map_free(Map* map) {
	if (map != NULL) {
		map_free(map->left);
		map_free(map->right);
		free(map->masks);
		free(map);
	}
}
/* map */


/* expression */
typedef struct ExprTree{
	value_t val;
	mask_t mask;
	struct ExprTree* left;
	op_t o;
	struct ExprTree* right;
} ExprTree;

ExprTree* exprtree_init(value_t val, mask_t mask) {
	ExprTree* expr = (ExprTree*)malloc(sizeof(ExprTree));
	if (expr == NULL) {
		fprintf(stderr, "Error: Memory allocation error (ExprTree)\n");
		exit(1);
	}

	expr->val = val;
	expr->mask = mask;
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
		exprtree_free(expr->left);
		exprtree_free(expr->right);
		free(expr);
	}
}
/* expression */


/* forest */
typedef struct ExprForest{
	ExprTree** expr;
	Map* map;
	mask_t maxmask;
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

	forest->map = map_init(0);
	forest->maxmask = 1;
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

		ExprTree** temp = (ExprTree**)realloc(forest->expr, forest->capacity * sizeof(ExprTree*));

		if (temp == NULL) {
			free(forest->expr);
			fprintf(stderr, "Error: Memory Reallocation error (exprforest_insert)\n");
			exit(1);
		}

		forest->expr = temp;
	}

	forest->expr[forest->top++] = expr;
}

void exprforest_input(ExprForest* forest, int N, int args[]) {
	int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};

	/* sort in desceding order of appearances */
	size_t u = 0;
	value_t nums[N];
	count_t copies[N];

	for (size_t i=0; i<N; ++i) {
		int place = -1;
		for (size_t j = 0; j<u; ++j)
			if (args[i] == nums[j])
				place = j;

		if (place == -1) {
			nums[u] = args[i];
			copies[u] = 1;
			++u;
		} else ++copies[place];
	}

	for (size_t i=0; i<u; ++i)
		for (size_t j=i+1; j<u; ++j)
			if (nums[i] < nums[j] && copies[i] < copies[j]) {
				value_t temp = nums[i]; nums[i] = nums[j]; nums[j] = temp;
				count_t temp2 = copies[i]; copies[i] = copies[j]; copies[j] = temp2;
			}

	for (size_t i=0; i<u; ++i) {
		ExprTree* expr = exprtree_init(nums[i], primes[i]);
		exprforest_insert(forest, expr);
		for (int j=0; j<copies[i]; ++j)
			forest->maxmask *= primes[i];
	}

	for (size_t i=0; i < forest->top; ++i)
		map_insert(forest->map, forest->expr[i]->val, forest->expr[i]->mask);

	forest->jumps = (size_t*)malloc((N+2) * sizeof(size_t));
	if (forest->jumps == NULL) {
		fprintf(stderr, "Error: Memory allocation error (jumps)\n");
		exit(1);
	}
	forest->jumps[1] = 0;
	forest->jumps[2] = u;
}

void exprforest_grow(ExprForest* forest, count_t target_num) {
	size_t top = forest->top;

	size_t next_break = 0;
	size_t next_break_idx = 2;
	count_t i_sum = 0;

	for (size_t i=0; i<top; ++i) {
		if (i >= next_break) {
			next_break = forest->jumps[next_break_idx++];
			++i_sum;
		}

		ExprTree* first = forest->expr[i];
		size_t L = forest->jumps[target_num - i_sum], U = forest->jumps[target_num - i_sum + 1];

		for (size_t j=L; j<U; ++j) {
			// if (first->mask > ULLONG_MAX / second->mask) continue;
			mask_t combmask = first->mask * forest->expr[j]->mask;
			if (forest->maxmask % combmask != 0)
				continue;

			ExprTree* second = forest->expr[j];

			/* Avoid repeated dereferencing */
			value_t first_val = first->val, second_val = second->val;
			op_t first_o = first->o, second_o = second->o;

			if (first_val < second_val || (first_val == second_val && first < second))
				continue;

			bool decreasing = first->right == NULL
			                || first->right->val > second -> val
			                || (first->right->val == second->val && first->right > second);

			op_t o = PLUS;
			value_t result = first_val + second_val;
			if (second_o != PLUS && second_o != MINUS
			   && ((first_o == TIMES || first_o == DIV || first_o == 0)
			       || (first_o == PLUS && decreasing))) {

				if (result < THRESHOLD_ && map_insert(forest->map, result, combmask)) {
					ExprTree* new_expr = exprtree_init(result, combmask);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					exprforest_insert(forest, new_expr);
				}
			}


			o = MINUS;
			result = first_val - second_val;
			if (result > 0 && second_o != PLUS && second_o != MINUS
			   && ((first_o == PLUS || first_o == TIMES || first_o == DIV || first_o == 0)
			       || (first_o == MINUS && decreasing))) {

				if (map_insert(forest->map, result, combmask)) {
					ExprTree* new_expr = exprtree_init(result, combmask);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					exprforest_insert(forest, new_expr);
				}
			}


			o = TIMES;
			result = first_val * second_val;
			if (second_o != TIMES && second_o != DIV
			   && ((first_o == PLUS || first_o == MINUS || first_o == 0)
			       || (first_o == TIMES && decreasing))) {

				if (result < THRESHOLD_ && map_insert(forest->map, result, combmask)) {
					ExprTree* new_expr = exprtree_init(result, combmask);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					exprforest_insert(forest, new_expr);
				}
			}


			o = DIV;
			result = first_val / second_val;
			if (first_val % second_val == 0 && second_o != TIMES && second_o != DIV
			   && ((first_o == PLUS || first_o == MINUS || first_o == TIMES || first_o == 0)
			       || (first_o == DIV && decreasing))) {

				if (map_insert(forest->map, result, combmask)) {
					ExprTree* new_expr = exprtree_init(result, combmask);
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
	forest->jumps[target_num + 1] = forest->top;
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
		printf("%llu\n", forest->expr[i]->mask);
		printf("\t");
		exprtree_printexpr(forest->expr[i]);
		printf("\n\n");
	}
}

void exprforest_free(ExprForest* forest) {
	for (size_t i=0; i < forest->top; ++i)
		free(forest->expr[i]);
	free(forest->expr);
	map_free(forest->map);
	free(forest->jumps);
	free(forest);
}
/* forest */


#endif /* EXPR_H_ */
