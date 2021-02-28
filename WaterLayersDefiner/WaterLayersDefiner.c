#include <Windows.h>
#include <stdio.h>

#pragma warning( disable : 4996 )

#define LINE_SIZE 256

int main(int argc, char * argv[])
{
	printf("Water layers definition tool is started\n");
	
	//Получение путей до файла PDB с полными данными моделирования и до файла, куда будут записываться очищенные данные

	char* path_to_mpf = argv[1];
	char* path_to_fpf = argv[2];

	printf("Path to main pdb data file: %s\nPath to filtered pdb data file: %s\n", path_to_mpf, path_to_fpf);

	FILE* main_pdb_data_file = fopen(path_to_mpf, "r");
	FILE* filtered_pdb_data_file = fopen(path_to_fpf, "r");

	char line[LINE_SIZE];
	int line_counter = 0;
	char ch = 0;
	while (ch = fgetc(main_pdb_data_file)) {
		printf("%c", ch);
		if (ch == '\n') {
			printf("\n");
			//break;
		}

	}
	printf("reading complete");


	//HANDLE* main_pdb_data_file = CreateFileA(path_to_mpf, GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL); 

	//HANDLE* filtered_pdb_data_file = CreateFileA(path_to_fpf, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);



	//Закрытие всех указателей на файлы
	//CloseHandle(main_pdb_data_file);
	//CloseHandle(filtered_pdb_data_file);
	
}