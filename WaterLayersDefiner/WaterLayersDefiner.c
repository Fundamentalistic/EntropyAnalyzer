#include <Windows.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#pragma warning( disable : 4996 )

#define LINE_SIZE 256

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

int ter_is_achieved;
int protein_size_val, protein_cursor;
atom* protein;

int protein_size(FILE* pdb);

int substring(char* str, int start, int finish, int linelen, char* result, int result_len)
{
	int result_cursor = 0;
	for (int i = start; i < finish; i++) {
		if (result_cursor >= result_len) {
			return 1;
		}
		result[result_cursor] = str[i];
		result_cursor++;
	}
	return 0;
}

int is_protein_record(char* line)
{
	char res[6], record[7];
	ZeroMemory(res, 6); ZeroMemory(record, 7);
	substring(line, 0, 6, LINE_SIZE, record, 6);
	substring(line, 17, 21, LINE_SIZE, res, 6);
	if (strstr(record, "ATOM") != 0){
		if (strstr(res, "TIP3") != 0) {
			return 0;
		}
		return 1;
	}
	return 0;
}

int is_end_of_frame(char* line)
{
	char record[7];
	substring(line, 0, 6, LINE_SIZE, record, 6);
	if (strstr(record, "END") != 0) {
		return 1;
	}
	return 0;
}

int protein_size(FILE* pdb) 
{
	char line[LINE_SIZE];
	int line_cursor = 0;
	char ch = 0;
	int protein_size_val = 0;
	int line_counter = 0;
	while (ch = fgetc(pdb)) {
		line[line_cursor] = ch;
		if (ch == '\n') {
			if( is_protein_record(line) ){
				protein_size_val++;
			}
			else if( is_end_of_frame(line) ) {
				return protein_size_val;
			}
			ZeroMemory(line, LINE_SIZE);
			line_counter++;
			line_cursor = 0;
			if (line_counter % 100000 == 0) {
				printf("%i lines was readed. Protein size: %i\r", line_counter, protein_size_val);
			}
			continue;
		}
		line_cursor++;
	}
}

int proceed_pdb_line(char* line, const int linelen)
{
	int
		id, resSeq, substring_status;

	double
		x, y, z;

	char
		idc[6], name[4], res[6], record[7], resSeqc[4], xc[8], yc[8], zc[8];

	ZeroMemory(idc, 6); ZeroMemory(name, 4); ZeroMemory(res, 6); ZeroMemory(record, 7);
	ZeroMemory(resSeqc, 4); ZeroMemory(xc, 4); ZeroMemory(yc, 5);  ZeroMemory(zc, 6);

	if (strstr(line, "END") != NULL) {
		ter_is_achieved = 1;
		return 2;
	}

	substring_status = substring(line, 7, 12, linelen, idc, 6);
	substring_status = substring(line, 0, 6, linelen, record, 6);
	substring_status = substring(line, 13, 17, linelen, name, 4);
	substring_status = substring(line, 17, 21, linelen, res, 6);
	substring_status = substring(line, 23, 27, linelen, resSeqc, 4);
	substring_status = substring(line, 31, 38, linelen, xc, 8);
	substring_status = substring(line, 38, 46, linelen, yc, 8);
	substring_status = substring(line, 46, 55, linelen, zc, 8);

	if (strstr(record, "ATOM") == NULL) {
		return 0;
	}
	id = atoi(idc);
	resSeq = atoi(resSeqc);
	x = atof(xc);
	y = atof(yc);
	z = atof(zc);
	name[3] = '\0';
	res[5] = '\0';
	if (strstr(res, "TIP3") == NULL) {
		if (!ter_is_achieved) {
			protein[protein_cursor].id = id;
			strcpy(protein[protein_cursor].name, name);
			strcpy(protein[protein_cursor].res, res);
			protein[protein_cursor].resSeq = resSeq;
			protein[protein_cursor].c.x = x;
			protein[protein_cursor].c.y = y;
			protein[protein_cursor].c.z = z;
			protein_cursor++;
		}
	}
	return 1;
}

int read_protein(FILE* pdb_file)
{
	char line[LINE_SIZE];
	int line_cursor = 0;
	char ch = 0;
	int protein_size_val = 0;
	int line_counter = 0;
	int res = 0;

	while (ch = fgetc(pdb_file)) {
		line[line_cursor] = ch;
		if (ch == '\n') {
			res = proceed_pdb_line(line, LINE_SIZE);
			if (res == 2) {
				return 0;
			}
			ZeroMemory(line, LINE_SIZE);
			line_counter++;
			line_cursor = 0;
			continue;
		}
		line_cursor++;
	}
}

int main(int argc, char * argv[])
{
	printf("Water layers definition tool is started\n");

	// Инициализация глобальных переменных

	ter_is_achieved = 0;
	protein_cursor = 0;
	
	//Получение путей до файла PDB с полными данными моделирования и до файла, куда будут записываться очищенные данные

	char* path_to_mpf = argv[1];
	char* path_to_fpf = argv[2];

	printf("Path to main pdb data file: %s\nPath to filtered pdb data file: %s\n", path_to_mpf, path_to_fpf);

	FILE* main_pdb_data_file = fopen(path_to_mpf, "r");
	if (main_pdb_data_file == NULL) {
		printf("Main pdb data file opening error: %i\n", GetLastError());
		return 1;
	}
	FILE* filtered_pdb_data_file = fopen(path_to_fpf, "w+");
	if (main_pdb_data_file == NULL) {
		printf("Filtered pdb data file opening error: %i\n", GetLastError());
		return 1;
	}

	// Получаем размер массива для хранения данных о белке

	protein_size_val = protein_size(main_pdb_data_file);

	// Выделяем память для белка

	protein = malloc(sizeof(atom) * protein_size_val);

	// Записываем операционные данные в структуру белка из файла PDB

	read_protein(main_pdb_data_file);
	printf("Protein cursor val: %i\n", protein_cursor);

	// Выполняем анализ 1, 2 и 3го водного слоя вокруг белка

	fclose(main_pdb_data_file);
	fclose(filtered_pdb_data_file);

	free(protein);
	
}