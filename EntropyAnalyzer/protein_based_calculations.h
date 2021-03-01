#ifndef PROTEIN_BASED
#define PROTEIN_BASED 1
#include <math.h>
#include "pdb.h"

#define M_PI 3.14159265358979323846

typedef struct {
	double sina;
	double cosa;
} step;

int atoms_count, atoms_cursor;
int water_size, water_cursor;
int protein_size, protein_cursor, mutex;
int active_threads;
#endif
