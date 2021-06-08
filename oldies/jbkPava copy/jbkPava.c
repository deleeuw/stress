#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

struct block {
	double value;
	double weight;
	int size;
	int previous;
	int next;
};

#define DEBUG false

void jbkPava (double *x, double *w, const int *n) {
	struct block *blocks = calloc ((size_t) * n, sizeof(struct block));
	for (int i = 0; i < *n; i++) {
		blocks[i].value = x[i];
		blocks[i].weight = w[i];
		blocks[i].size = 1;
		blocks[i].previous = i - 1; // index first element previous block
		blocks[i].next = i + 1;     // index first element next block
	}
	int active = 0;
	do {
		bool upsatisfied = false;
		int next = blocks[active].next; 
		if (next == *n) upsatisfied = true;
		else if (blocks[next].value > blocks[active].value) upsatisfied = true;
		if (!upsatisfied) {
			double ww = blocks[active].weight + blocks[next].weight;
			int nextnext = blocks[next].next;
			blocks[active].value = (blocks[active].weight * blocks[active].value + blocks[next].weight * blocks[next].value) / ww;
			blocks[active].weight = ww;
			blocks[active].size += blocks[next].size;
			blocks[active].next = nextnext;
			if (nextnext < *n)
				blocks[nextnext].previous = active;
			blocks[next].size = 0;
		}
		if (DEBUG) {
			printf("active = %2d upsatisfied = %2d\n", active, upsatisfied);
			for (int i = 0; i < *n; i++) {
				printf ("%8.4f %8.4f %2d %2d %2d\n", blocks[i].value, blocks[i].weight, blocks[i].size, blocks[i].previous, blocks[i].next);
			}
			printf("+++++++++++++++++++++++++++++++\n");
		}
		bool downsatisfied = false;
		int previous = blocks[active].previous;
		if (previous == -1) downsatisfied = true;
		else if (blocks[previous].value < blocks[active].value) downsatisfied = true;
		if (!downsatisfied) {
			double ww = blocks[active].weight + blocks[previous].weight;
			int previousprevious = blocks[previous].previous;
			blocks[active].value = (blocks[active].weight * blocks[active].value + blocks[previous].weight * blocks[previous].value) / ww;
			blocks[active].weight = ww;
			blocks[active].size += blocks[previous].size;
			blocks[active].previous = previousprevious;
			if (previousprevious > -1)
				blocks[previousprevious].next = active;
			blocks[previous].size = 0;
		}
		if (DEBUG) {
			printf("active = %2d downsatisfied = %2d\n", active, downsatisfied);
			for (int i = 0; i < *n; i++) {
				printf ("%8.4f %8.4f %2d %2d %2d\n", blocks[i].value, blocks[i].weight, blocks[i].size, blocks[i].previous, blocks[i].next);
			}
			printf("+++++++++++++++++++++++++++++++\n");
		}
		if ((blocks[active].next == *n) && downsatisfied) break;
		if (upsatisfied && downsatisfied) active = next;
	} while (true);
	int k = 0;
	for (int i = 0; i < *n; i++) {
		int blksize = blocks[i].size;
		if (blksize > 0.0) {
			for (int j = 0; j < blksize; j++) {
				if (DEBUG) printf("%2d %2d %2d %2d\n", i, j, k, blksize);
				x[k] = blocks[i].value;
				k++;
			}
		}
	}
	free (blocks);
}