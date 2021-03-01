#ifndef PDB
#define PDB 1
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include "protein_based_calculations.h"

#pragma warning( disable : 4996 )

typedef struct {
	double x;
	double y;
	double z;
} vector;

typedef struct {
	double p;
	double phi;
	double z;
} cylindric_point;

typedef struct {
	int id;
	char name[4];
	char res[6];
	int resSeq;
	vector c;
	cylindric_point cp;
	char for_calculate;
} atom;

atom* protein;
atom* water;

int read_pdb_file(char* path);
int proceed_pdb_line(char* line, const int linelen);
int str_is(char* str1, int str1len, const char* str2, int str2len);
int prepare_memory_for_data_storage(void);
int substring(char* str, int start, int finish, int linelen, char* result, int result_len);
int calculate_water_orientational_distribution(vector* CA, vector* CB, vector* N, atom* HOHa, int HOHa_len, double dist, FILE* log);
void rotation_x(vector* v, double cosa, double sina);
void rotation_y(vector* v, double cosa, double sina);
void rotation_z(vector* v, double cosa, double sina);
void shift(vector* shift, vector* vector);
void free_pdb_memory_stack();
void print_vector(vector* a, const char* description);
void proceed_pdb_data(void);
double distance(vector* a, vector* b);
double check_if_null(double val);
double cos_a(vector* a, vector* b);
double sin_from_cos(double cosa);

#endif
