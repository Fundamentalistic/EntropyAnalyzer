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
	if (strstr("ATOM", record) != 0) {
		if (strstr("TIP3W", res) != 0) {
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
	if (strstr("END", record) != 0) {
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
			line_cursor = 0;
		}
	}
}

int main(int argc, char * argv[])
{
	printf("Water layers definition tool is started\n");
	
	//Получение путей до файла PDB с полными данными моделирования и до файла, куда будут записываться очищенные данные

	char* path_to_mpf = argv[1];
	char* path_to_fpf = argv[2];

	printf("Path to main pdb data file: %s\nPath to filtered pdb data file: %s\n", path_to_mpf, path_to_fpf);

	FILE* main_pdb_data_file = fopen(path_to_mpf, "r");
	FILE* filtered_pdb_data_file = fopen(path_to_fpf, "r");

	int protein_size_val = protein_size(main_pdb_data_file);
	printf("Protein size: %i\n", protein_size_val);


	//HANDLE* main_pdb_data_file = CreateFileA(path_to_mpf, GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL); 

	//HANDLE* filtered_pdb_data_file = CreateFileA(path_to_fpf, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);



	//Закрытие всех указателей на файлы
	//CloseHandle(main_pdb_data_file);
	//CloseHandle(filtered_pdb_data_file);
	
}