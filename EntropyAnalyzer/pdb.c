#include "pdb.h"
#include <stdio.h>

#define PDB_LINE_SIZE 79

char* reading_buffer = NULL;
DWORD64 file_size;
atoms_count = 0, atoms_cursor = 0;
water_size = 0, water_cursor = 0;
protein_size = 0, protein_cursor = 0;
int line_cursor = 0;
int ter_is_achieved = 0;



int read_pdb_file(char* path)
{
	printf("PDB %s READING ...\n", path);
	HANDLE* file = CreateFileA(path, GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (file == INVALID_HANDLE_VALUE) {
		printf("PDB OPENING ERROR: %i\n", GetLastError());
		return 0;
	}
	BOOL res = GetFileSizeEx(file, &file_size); 
	printf("FILESIZE: %lli\n", file_size);
	if (file_size == -1) {
		printf("PDB FILE SIZE\n");
		return 0;
	}
	long buffer_cursor = 0;
	reading_buffer = (char*)malloc(file_size);
	DWORD64 readed = 0;
	res = ReadFile(file, reading_buffer, file_size, &readed, NULL);
	if (res == 0) {
		printf("Reading pdb file error: %i\n", GetLastError());
		return 0;
	}
	printf("READED: %lli\n", readed);
	CloseHandle(file);
	return 1;
}

char line[128];

void proceed_pdb_data(void)
{
	ter_is_achieved = 0;
	int line_cursor = 0;
	ZeroMemory(line, 128);
	double percents = 0;
	for (DWORD64 i = 0; i < file_size; i++) {
		if (i % 1000000 == 0) {
			percents = (double)i * 100 / file_size;
			printf("\tproceed_pdb_data: %lf\r", percents);
		}
		if (reading_buffer[i] == '\n') {
			proceed_pdb_line(line, 128);
			ZeroMemory(line, 128);
			line_cursor = 0;
		}
		else {
			line[line_cursor] = reading_buffer[i];
			line_cursor++;
			if (line_cursor >= 128) {
				line_cursor = 0;
			}
		}
	}
	free(reading_buffer);
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
		return 1;
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
	if (strstr(res, "TIP3") != NULL) {
		water[water_cursor].id = id;
		strcpy(water[water_cursor].name, name);
		strcpy(water[water_cursor].res, res);
		water[water_cursor].resSeq = resSeq;
		water[water_cursor].c.x = x;
		water[water_cursor].c.y = y;
		water[water_cursor].c.z = z;
		water_cursor++;
	}
	else {
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

int str_is(char* str1, int str1len, const char* str2, int str2len)
{
	if (str1len == str2len) {
		for (int i = 0; i < str1len; i++) {
			if (str1[i] != str2[i]) {
				return 0;
			}
		}
	}
	else {
		return 0;
	}
	return 1;
}

int prepare_memory_for_data_storage(void)
{
	char line[128];
	int line_cursor = 0;
	ZeroMemory(line, 128);
	ter_is_achieved = 0;
	double percents = 0;
	for (DWORD64 i = 0; i < file_size; i++) {
		if (i % 1000000 == 0) {
			percents = (double)(i * 100 / file_size);
			printf("\tprepare_memory_for_data_storage: %lf\r", percents);
		}
		if (reading_buffer[i] == '\n') {
			if (strstr(line, "ATOM  ") != NULL) {
				atoms_count++;
				if (strstr(line, "TIP3") != NULL) {
					water_size++;
				}
				else {
					if (!ter_is_achieved) {
						protein_size++;
					}
				}
			}
			if (strstr(line, "END") != NULL) {
				ter_is_achieved = 1;
			}
			ZeroMemory(line, 128);
			line_cursor = 0;
		}
		else {
			line[line_cursor] = reading_buffer[i];
			line_cursor++;
			if (line_cursor >= 128) {
				line_cursor = 0;
			}
		}
	}
	//printf("\r                                                                              \r");
	printf("protein_size: %i, water_size: %i\n", protein_size, water_size);
	protein = (atom*)malloc(sizeof(atom) * protein_size);
	water = (atom*)malloc(sizeof(atom) * water_size);

	return 1;
}

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

void free_pdb_memory_stack()
{
	free(water);
	free(protein);
}





