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
#include<string.h>
#include<stdbool.h>
#include<omp.h>

#include<limits.h>
#include<sys/resource.h>

#define START_SIZE_FOREST_ 200
#define START_SIZE_MAT_ 50

#define THRESHOLD_ ULLONG_MAX

#define MAX_MEMORY_ 9*1024*1024  /* in KiB */

#define p 196613 /* hashing prime */

typedef enum op_t {
	UNDEF = 0,
	PLUS = 43,
	MINUS = 45,
	TIMES = 42,
	DIV = 47,
} op_t;

typedef enum color_t {
	black, red
} color_t;

typedef unsigned long long value_t;
typedef unsigned char count_t;
typedef unsigned long long mask_t;

/* struct ExprTree 
 * represents arithmetic expression as binary tree
 * val ... final value
 * left, right pointers to left/right operand
 * o ... operation between operands


 * ExprTree* exprtree_init(value_t, mask_t) ... initializer
 * void exprtree_printexpr(const ExprTree*) ... print in infix form
 * void exprtree_free(ExprTree*) ... free pointers
*/
typedef struct ExprTree{
	value_t val;
	struct ExprTree* left;
	op_t o;
	struct ExprTree* right;
} ExprTree;

ExprTree* exprtree_init(value_t val) {
	ExprTree* expr = (ExprTree*)malloc(sizeof(ExprTree));
	if (expr == NULL) {
		fprintf(stderr, "Error: Memory allocation error (ExprTree)\n");
		exit(1);
	}

	expr->val = val;
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
/* ExprTree */

/* Cluster - array of ExprTrees with the same mask
 * size, capacity ... of container
 * exprs ... array of ExprTree*
 *
 * Cluster* cluster_init(mask_t, size_t) ... initialize with mask and expected maximum size
 * void cluster_insert(Cluster*, ExprTree*) ... add ExprTree* to cluster
 * void cluster_free(Cluster*) - free memory
*/
typedef struct Cluster {
	size_t size, capacity;
	mask_t mask;
	ExprTree** exprs;
} Cluster;

Cluster* cluster_init(mask_t mask, size_t expected_size) {
	Cluster* cluster = (Cluster*)malloc(sizeof(Cluster));
	if (cluster == NULL) {
		fprintf(stderr, "Malloc error in cluster_init.\n");
		exit(1);
	}
	
	ExprTree** exprs = (ExprTree**)malloc(expected_size * sizeof(ExprTree*));
	if (exprs == NULL) {
		fprintf(stderr, "Malloc error in cluster_init.\n");
		exit(1);
	}
	
	cluster->size = 0;
	cluster->capacity = expected_size;
	cluster->mask = mask;
	cluster->exprs = exprs;
	return cluster;
}

void cluster_insert(Cluster* cluster, ExprTree* expr) {
	if (cluster->size == cluster->capacity) {
		cluster->capacity *= 2;

		ExprTree** temp = (ExprTree**)realloc(cluster->exprs, cluster->capacity * sizeof(ExprTree*));
		
		if (temp == NULL) {
			free(cluster->exprs);
			fprintf(stderr, "Realloc error in cluster_insert.\n");
			exit(1);
		}

		cluster->exprs = temp;
	}

	cluster->exprs[cluster->size++] = expr;
}

void cluster_free(Cluster* cluster) {
	for (size_t i=0; i < cluster->size; ++i)
		free(cluster->exprs[i]);
	free(cluster->exprs);
	free(cluster);
}
/* Cluster */


/* Map - red-and-black search tree with expression results (val) as data
 * to each value we map an array of masks
 * Map* left, right, parent
 * color_t color
 *
 * value_t val 
 * mask_t* masks ... the masks with which val can be achieved
 * size, capacity ... of container
 *
 * Map* map_init(value_t, color_t) ... initialize new tree with given value and color
 * void map_left_rotate(Map*) ... left rotate on red-and-black tree
 * void map_right_rotate(Map*) ... right rotate on red-and-black tree
 * void map_rebalance(Map*) ... rebalancing red-and-black tree after insertion
 * bool map_insert_mask(Map*, mask_t) ... insert new mask if there is no divisor, return 0 if suboptimal
 * bool map_insert(Map*, value_t, mask_t) ... try to insert new value with mask, return 1 if inserted
 * void map_free(Map*) ... free memory
*/
typedef struct Map{
	struct Map *left, *right, *parent;
	color_t color;

	value_t val;

	mask_t* masks;
	size_t size;
	size_t capacity;
} Map;



Map* map_init(value_t val, color_t color) {
	Map* map = (Map*)malloc(sizeof(Map));
	if (map == NULL) {
		fprintf(stderr, "Error: Memory allocation error (map_init)\n");
		exit(1);
	}

	map->left = map->right = map->parent = NULL;
	map->color = color;
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

void map_left_rotate(Map* x) {
	if (x == NULL || x->right == NULL)
		return;
	
	Map* y = x->right;

	x->right = y->left;
	if (x->right != NULL) 
		x->right->parent = x;
	
	y->parent = x->parent;
	if (x->parent != NULL)
		if (x == x->parent->left)
			x->parent->left = y;
		else
			x->parent->right = y;

	y->left = x;
	x->parent = y;
}

void map_right_rotate(Map* y) {
	if (y == NULL || y->left == NULL)
		return;

	Map* x = y->left;

	y->left = x->right;
	if (y->left != NULL)
		y->left->parent = y;
	
	x->parent = y->parent;
	if (y->parent != NULL)
		if (y == y->parent->left)
			y->parent->left = x;
		else
			y->parent->right = x;

	x->right = y;
	y->parent = x;
}

void map_rebalance(Map* map) {
	Map *parent = NULL, *grandparent = NULL, *uncle = NULL;

	while (map->parent != NULL && map->color == red && map->parent->color == red) {
		parent = map->parent;
		grandparent = parent->parent;

		if (parent == grandparent->left)
			uncle = grandparent->right;
		else
			uncle = grandparent->left;

		if (uncle == NULL)
			map = grandparent;
		else if (uncle->color == red) {
			uncle->color = black;
			parent->color = black;
			grandparent->color = red;
			map = grandparent;
		}

		else {
			/* LL */
			if (parent == grandparent->left && map == parent->left) {
				color_t temp = parent->color;
				parent->color = grandparent->color;
				grandparent->color = temp;

				map_right_rotate(grandparent);
			}
			
			/* LR */
			else if (parent == grandparent->left && map == parent->right) {
				color_t temp = map->color;
				map->color = grandparent->color;
				grandparent->color = temp;

				map_left_rotate(parent);
				map_right_rotate(grandparent);
			}

			/* RR */
			else if (parent == grandparent->right && map == parent->right) {
				color_t temp = parent->color;
				parent->color = grandparent->color;
				grandparent->color = temp;

				map_left_rotate(grandparent);
			}

			/* RL */
			else {
				color_t temp = map->color;
				map->color = grandparent->color;
				grandparent->color = temp;

				map_right_rotate(parent);
				map_left_rotate(grandparent);
			}
		}
	}		
}

static inline bool map_insert_mask(Map* map, mask_t mask) {
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

static inline bool map_insert(Map* map, value_t val, mask_t mask) {
	Map* new_map;
	int count = 0;

	while (1) {
		++count;
		if (map->val == val) { //printf("I was here .... %d\n", count);
			return map_insert_mask(map, mask);} 

		else if (val < map->val) {
			if (map->left != NULL)
				map = map->left;
			else {
				new_map = map_init(val, red);
				map->left = new_map;
				new_map->parent = map;
				map_rebalance(new_map);
//				map_insert_mask(new_map, mask);	
			}
		}
		else {
			if (map->right != NULL) 
				map = map->right;
			else {
				new_map = map_init(val, red);
				map->right = new_map;
				new_map->parent = map;
				map_rebalance(new_map);
//				map_insert_mask(new_map, mask); //return?
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
		printf("Map at address %p has color %d and children...\n", map, map->color);
		printf("%llu has masks:\n", map->val);
		for (int i=0; i < map->size; ++i)
			printf("%llu ", map->masks[i]);
		printf("\n");
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

/* HashTable - simple mod based hashtable
 * maps masks to indices
 *
 * mask_t[] masks ... masks column
 * int[] indices ... indices column (-1 if none is mapped)
 *
 * HashTable* ht_init() ... initialize new hashtable
 * void ht_save_idx(HashTable*, mask_t, size_t) ... add mask and index
 * int ht_get_idx(HashTable*, mask_t) ... get index of mask or -1
 *
*/

typedef struct HashTable {
	mask_t masks[p];
	int indices[p];
} HashTable;

HashTable* ht_init() {
	HashTable* ht = (HashTable*)malloc(sizeof(HashTable));

	if (ht == NULL) {
		fprintf(stderr, "Memory allocation error (ht_init).\n");
		exit(1);
	}

	memset(ht->masks, 0, p * sizeof(mask_t));
	memset(ht->indices, -1, p * sizeof(int));
	return ht;
}

void ht_save_idx(HashTable* ht, mask_t mask, size_t idx) {
	size_t loc = mask % p;
	while (ht->masks[loc]) 
		if (++loc == p) loc = 0;

	ht->masks[loc] = mask;
	ht->indices[loc] = idx;
}

int ht_get_idx(HashTable* ht, mask_t mask) {
	size_t loc = mask % p;
	while (ht->masks[loc] != mask && ht->masks[loc] != 0)
		if (++loc == p) loc = 0;
	
	return ht->indices[loc];
}

/* HashTable*/

/* ExprForest ... main structure collecting and combining expressions
 * 
 * Cluster** clusters ... array of clusters
 * HashTable* ht ... hashtable mapping masks to indices in clusters
 * Map* map ... maps masks to achieved results
 * mask_t maxmask ... represents the whole set of numbers; any other mask must divide it
 * size_t* jumps ... array representing boundaries of cluster groups (each group uses same set size)
 * size_t top, capacity ... of clusters' container
 *
 * ExprForest* exprforest_init() ... prepare empty forest, allocate memory
 * void exprforest_insert(ExprForest*, Cluster*) ... add cluster to array
 * void exprforest_input(ExprForest*, int, int[])
 * 	... remove duplicates, sort values, assign base masks
 * 	... initialize maxmask, map, hashtable, clusters and jumps
 *
 * void exprforest_grow(ExprForest*, count_t) ... combine clusters to get all clusters of target size
 *
 * int exprforest_count(const ExprForest*) ... return total number of solutions
 * void exprforest_query(const ExprForest*, value_t, bool) ... print all expressions with target value
 * 
 * void exprforest_free(ExprForest*) ... deallocate memory
*/
typedef struct ExprForest{
	Cluster** clusters;
	HashTable* ht;
	Map* map;
	mask_t maxmask;
	size_t* jumps;
	size_t top, capacity;
} ExprForest;

ExprForest* exprforest_init() {
	ExprForest* forest = (ExprForest*)malloc(sizeof(ExprForest));
	if (forest == NULL) {
		fprintf(stderr, "Error: Memory allocation error (exprforest)\n");
		exit(1);
	}

	forest->clusters = (Cluster**)malloc(START_SIZE_FOREST_ * sizeof(Cluster*));
	if (forest->clusters == NULL) {
		fprintf(stderr, "Error: Memory allocation error(exprforest expr)\n");
		exit(1);
	}
	
	forest->ht = ht_init();

	forest->map = map_init(0, black);
	forest->maxmask = 1;
	forest->jumps = NULL;
	forest->top = 0;
	forest->capacity = START_SIZE_FOREST_;

	return forest;
}

void exprforest_insert(ExprForest* forest, Cluster* cluster) {
	if (forest->top == forest->capacity) {
		forest->capacity *= 2;

		Cluster** temp = (Cluster**)realloc(forest->clusters, forest->capacity * sizeof(Cluster*));

		if (temp == NULL) {
			free(forest->clusters);
			fprintf(stderr, "Error: Memory Reallocation error (exprforest_insert)\n");
			exit(1);
		}

		forest->clusters = temp;
	}

	forest->clusters[forest->top++] = cluster;
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

	
	/* initialize starting clusters, expressions, hashtable and map */
	for (size_t i=0; i<u; ++i) {
		Cluster* cluster = cluster_init(primes[i], 1);
		ExprTree* expr = exprtree_init(nums[i]);

		cluster_insert(cluster, expr);
		exprforest_insert(forest, cluster);

		ht_save_idx(forest->ht, primes[i], i);

		map_insert(forest->map, nums[i], primes[i]);

		for (int j=0; j<copies[i]; ++j)
			forest->maxmask *= primes[i];
	}
	
	/* initialize jumps */
	forest->jumps = (size_t*)malloc((N+1) * sizeof(size_t));
	if (forest->jumps == NULL) {
		fprintf(stderr, "Error: Memory allocation error (jumps)\n");
		exit(1);
	}

	/* jumps[i] is starting index for clusters of size i+1 */
	forest->jumps[0] = 0;
	forest->jumps[1] = u;
}

void exprforest_free(ExprForest*);

void exprforest_grow(ExprForest* forest, count_t target_num) {
	/* freeze top */
	size_t top = forest->top;
	
	/* I_sum keeps track of current cluster size */
	size_t next_break = 0;
	size_t next_break_idx = 1;
	count_t I_sum = 0;

	int target_idx = -1;

	Cluster *first_cluster, *second_cluster, *target_cluster;
	ExprTree *first, *second;
	
	/* iterate among all clusters
	 * with target_num for combined size*/
	for (size_t I=0; I<top; ++I) {
		if (I >= next_break) {
			next_break = forest->jumps[next_break_idx++];
			++I_sum;
		}
		
		first_cluster = forest->clusters[I];
		size_t L = forest->jumps[target_num - I_sum - 1], U = forest->jumps[target_num - I_sum];

		for (size_t J=L; J<U; ++J) {
			second_cluster = forest->clusters[J];
			mask_t combmask = first_cluster->mask * second_cluster->mask;
			if (forest->maxmask % combmask != 0)
				continue;
			
			/* insert or set target cluster */
			if ((target_idx = ht_get_idx(forest->ht, combmask)) == -1) {
				target_cluster = cluster_init(combmask, first_cluster->size * second_cluster->size * 2);
				ht_save_idx(forest->ht, combmask, forest->top);
				exprforest_insert(forest, target_cluster);
			} else target_cluster = forest->clusters[target_idx];
			
			value_t result, first_val, second_val;
			op_t o, first_o, second_o;
			/* try to combine all pairs of expressions */
			for (size_t i=0; i < first_cluster->size; ++i)
			for (size_t j=0; j < second_cluster->size; ++j) {
			first = first_cluster->exprs[i];
			second = second_cluster->exprs[j];

			/* avoid repeated dereferencing */
			first_val = first->val, second_val = second->val;
			first_o = first->o, second_o = second->o;
		
			/* skip if first expr is smaller; the same pair will come again in reversed order */
			if (first_val < second_val || (first_val == second_val && first < second))
				continue;

			/* in case of chained operations arguments must be increasing */
			bool decreasing = first->right == NULL
			                || first->right->val > second -> val
			                || (first->right->val == second->val && first->right > second);
			
			/* for each of +, -, *, / check conditions
			 * if they hold, insert into target_cluster */
			o = PLUS;
			result = first_val + second_val;
			if (result < THRESHOLD_ && second_o != PLUS && second_o != MINUS
			   && ((first_o == TIMES || first_o == DIV || first_o == 0)
			       || (first_o == PLUS && decreasing))) {

				if (map_insert(forest->map, result, combmask)) {
					ExprTree* new_expr = exprtree_init(result);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					cluster_insert(target_cluster, new_expr);
				}
			}

			o = MINUS;
			result = first_val - second_val;
			if (result > 0 && second_o != PLUS && second_o != MINUS
			   && ((first_o == PLUS || first_o == TIMES || first_o == DIV || first_o == 0)
			       || (first_o == MINUS && decreasing))) {
				
				if (map_insert(forest->map, result, combmask)) {
					ExprTree* new_expr = exprtree_init(result);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					cluster_insert(target_cluster, new_expr);
				}
			}

			o = TIMES;
			result = first_val * second_val;
			if (result < THRESHOLD_ && second_o != TIMES && second_o != DIV
			   && ((first_o == PLUS || first_o == MINUS || first_o == 0)
			       || (first_o == TIMES && decreasing))) {

				if (map_insert(forest->map, result, combmask)) {
					ExprTree* new_expr = exprtree_init(result);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					cluster_insert(target_cluster, new_expr);
				}
			}

			o = DIV;
			result = first_val / second_val;
			if (first_val % second_val == 0 && second_o != TIMES && second_o != DIV
			   && ((first_o == PLUS || first_o == MINUS || first_o == TIMES || first_o == 0)
			       || (first_o == DIV && decreasing))) {

				if (map_insert(forest->map, result, combmask)) {
					ExprTree* new_expr = exprtree_init(result);
					new_expr->left = first;
					new_expr->right = second;
					new_expr->o = o;
					cluster_insert(target_cluster, new_expr);
				}
			}
		}}

		if (I%500 == 0) {
			struct rusage usage;
			getrusage(RUSAGE_SELF, &usage);
			if (usage.ru_maxrss > MAX_MEMORY_) {
				fprintf(stderr, "Memory limit exceeded (%d MiB). The program will terminate.\n", MAX_MEMORY_/1024);
				exprforest_free(forest);
				exit(1);
			}
		}
	}
	forest->jumps[target_num] = forest->top;

	printf("I have %d clusters. Coming with masks and sizes ...\n", forest->top);
	for (int i=0; i<forest->top; ++i) {
		printf("Index %d Mask %d (quot %d )Size %d\n",
				i, forest->clusters[i]->mask,
				forest->maxmask/forest->clusters[i]->mask,
				forest->clusters[i]->size);
	}
}


int exprforest_count(const ExprForest* forest) {
	int count = 0;
	for (size_t i=0; i < forest->top; ++i)
		count += forest->clusters[i]->size;

	return count;
}

void exprforest_query(const ExprForest* forest, value_t query) {
	int count = 0;
	printf("\nQuery: %llu\n", query);
	for (size_t i=0; i < forest->top; ++i)
		for (size_t j=0; j < forest->clusters[i]->size; ++j)
			if (forest->clusters[i]->exprs[j]->val == query) {
				exprtree_printexpr(forest->clusters[i]->exprs[j]);
				printf("\n");
				++count;
			}

	printf("%d solutions found for %llu.\n------\n", count, query);
}

void exprforest_free(ExprForest* forest) {
	for (size_t i=0; i < forest->top; ++i)
		cluster_free(forest->clusters[i]);
	free(forest->clusters);
	free(forest->ht);
	map_free(forest->map);
	free(forest->jumps);
	free(forest);
}
/* forest */


#endif /* EXPR_H_ */
