/*
* Анализатор водной ориентации вокруг полпептидной цепи
* Автор: Беспалов Андрей Витальевич. Аспирант ФБУН ГНЦ ВБ "Вектор"
* Программа принимает на вход путь к pdb файлу содержащему кадры моделирования молекулярной динамики
* Выход программы составляет ориентационное распределение воды соответствующее каждой аминокислоте системы
*/

/*
* Программа не рассчитана на чтение файлов, объем которых превышает объем оперативной памяти
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <windows.h>
#include "constants.h"
#include "pdb.h"
#include "protein_based_calculations.h"  // Содержит код функции calculate_water_orientational_distribution,
// которая вычисляет ориентационноее распределение и пишет файл log.

#define CA_ATOM 0
#define CB_ATOM 1
#define N_ATOM 2
#define C_ATOM 3
#define OH2_ATOM 4
#define H1_ATOM 5
#define H2_ATOM 6
#define TIP3_RES 1
#define NON_WATER_RES 2

#define CA 0
#define CB 1
#define N 2
#define C 3

#define LOGPATH "C:\\EntropyLogs\\log%i.txt"
#define LOG_FILES_STRING "C:\\EntropyLogs\\log%i.txt"

#define DISTANCE 50

void reset_for_calculate();
void generate_cylindrical_coordinates(atom* n, atom* c);

int main(int argc, char* argv[]) {

	if (argc < 2) {
		printf("FIRST_ARGUMENT: %s", argv[0]);
		printf(TUTORIAL);
		return 0;
	}

	char* link_on_path_argument = argv[1];
	int res = 0;
	printf("try to read: %s\n", link_on_path_argument);
	res = read_pdb_file(link_on_path_argument);
	if (res == 0) {
		printf("Reading file error");
		return 1;
	}
	printf("try to prepare memory: %s\n", link_on_path_argument);
	prepare_memory_for_data_storage();
	printf("proceed pdb data : %s\n", link_on_path_argument);
	proceed_pdb_data();

	atom* calculation_sequence = (atom*)malloc(sizeof(atom) * 4);

	ZeroMemory(calculation_sequence, (sizeof(atom) * 4));
	calculation_sequence[0].id = -1;
	calculation_sequence[1].id = -1;
	calculation_sequence[2].id = -1;
	calculation_sequence[3].id = -1;

	int calculation_result = 0, log_index = 0, percents = 0;
	char logname[512];

	for (int i = 0; i < protein_size; i++) { // для всех атомов белка

		if (strstr("CA ", protein[i].name) != 0) {
			// если имя атома СА
			calculation_sequence[CA] = protein[i];
		}
		else if (strstr("CB ", protein[i].name) != 0 || strstr("HA1", protein[i].name) != 0) {
			// если имя атома СВ или НА (боковой радикал глицина)
			calculation_sequence[CB] = protein[i];
		}
		else if (strstr("N  ", protein[i].name) != 0) {
			// если имя атома N
			calculation_sequence[N] = protein[i];
		}
		else if (strstr("C  ", protein[i].name) != 0) {
			// если имя атома C
			calculation_sequence[C] = protein[i];
		}
		if (
			calculation_sequence[CA].id != -1 &&
			calculation_sequence[CB].id != -1 &&
			calculation_sequence[N].id != -1 &&
			calculation_sequence[C].id != -1
			) {
			// если все атомы вычисляемой последовательности заданы, начинаем вычисление ориентационного распределения
			printf("Generate cylindrical coordinates for each point\n");
			generate_cylindrical_coordinates(&calculation_sequence[N_ATOM], &calculation_sequence[C_ATOM]);
			ZeroMemory(logname, 512);
			printf(LOG_FILES_STRING, log_index);
			sprintf(logname, LOG_FILES_STRING, log_index);
			FILE* log_file = fopen(logname, "a+");

			if (!log_file) {
				printf("%s opening error\n", logname);
				return 1;
			}
			else {
				printf("Fopen result: %i, %s\n", log_file, logname);
			}
			// Требуется метод выборки массива воды из диапазона аппликат
			percents = i * 100 / protein_size;
			printf("COMPLETE: %i, LOGFILE: %s\n", percents, logname);

			calculation_result = calculate_water_orientational_distribution(
				&calculation_sequence[CA].c, // Координаты СА атома
				&calculation_sequence[CB].c, // Координаты СВ атома
				&calculation_sequence[N].c,  // Координаты N атома
				water,						 // Массив молекул воды
				water_size,					 // Размер массива молекул воды
				DISTANCE,					 // Дистанция на которой учитывается ориентационное распределение
				log_file				     // Дескриптор файла лога
			);
			log_index++;
			reset_for_calculate();
			if (log_file) {
				fclose(log_file);
			}
			if (calculation_result == -1) {
				printf("CALCULATION ERROR. PROCEDURE INTERRUPTION");
				break;
			}

			calculation_sequence[0].id = -1;
			calculation_sequence[1].id = -1;
			calculation_sequence[2].id = -1;
			calculation_sequence[3].id = -1;
		}
	}

	/*
		1) Пробегаем про файлу PDB, индекс лога равен 0

		Для вычисления каждого файла логов, который в последствии станет объектом DataFrame, нам необходимо в начале вычисления иметь следующие данные
		1.1 координаты С, СА, СВ, N атомов
		1.2 массив молекул воды, его размер

		2) Находим направляющий вектор прямой, проходящей через C и N атомы
		2.1 Совмещаем ось Z с найденным на 2 направляющим вектором и N с началом декартовой	системы координат
		2.2 Находим угол между проекцей CA на XOY и осью X
		2.3 Выполняем поворот вокруг Z на найденный угол
		2.4 Преобразовываем координаты всех атомов вычисляемой системы к новым координатам после поворота и сдвига
		2.5 Преобразовываем все координаты атомов вычисляемой системы к цилиндрическим координатам
		3) находим все атомы кислорода между Nz и Cz
		3.1) Выполняем копирование указателей на участвующие в вычислении атомы
		3.2) Выполняем вычисление ориентационного распределения воды для данного массива молекул воды
		3.3) Записываем результат в файл лога с индексом лога
		3.4) Индекс лога увеличиваем на единицу
	*/

	free_pdb_memory_stack();
	free(calculation_sequence);
	return 0;
}

void reset_for_calculate()
{
	for (int i = 0; i < water_size; i++) {
		water[i].for_calculate = 0;
	}
}

void generate_cylindrical_coordinates(atom* n, atom* c)
{
	/*
	* Переводим декартовы координаты в новую систему координат с совмещенными осями Z и направляющим вектором
	*/

	// Вычисляем направляющий вектор
	vector z_dir;
	z_dir.x = c->c.x - n->c.x;
	z_dir.y = c->c.y - n->c.y;
	z_dir.z = c->c.z - n->c.z;

	// Вычисляем косинус угла между проекцией направляющего ветора на ХОZ

	step steps[2];

	vector z_dir_xoz_proj;
	z_dir_xoz_proj.x = z_dir.x;
	z_dir_xoz_proj.y = 0;
	z_dir_xoz_proj.z = z_dir.z;

	vector z;
	z.x = 0; z.y = 0; z.z = 1;

	double cosa = cos_a(&z_dir_xoz_proj, &z);
	double alpha = 0;
	double hemi_pi = M_PI / 2;

	if (z_dir_xoz_proj.z < 0) {
		alpha = acos(cosa);
		cosa = cos((hemi_pi - alpha) + hemi_pi) * -1;
	}
	double sina = sin_from_cos(cosa);
	if (z_dir_xoz_proj.x > 0) {
		sina = -sina;
	}
	steps[0].sina = sina;
	steps[0].cosa = cosa;
	rotation_y(&z_dir, cosa, sina);
	cosa = cos_a(&z_dir, &z);
	if (z_dir.x < 0) {
		printf("calculation error. negative Z");
	}
	sina = sin_from_cos(cosa);
	if (z_dir.y < 0) {
		sina = -sina;
	}

	steps[1].sina = sina;
	steps[1].cosa = cosa;

	rotation_x(&z_dir, cosa, sina);

	vector h;

	for (int i = 0; i < water_size; i++) {
		if (strcmp(water[i].name, "OH2") == 0) {
			h.x = water[i].c.x;
			h.y = water[i].c.y;
			h.z = water[i].c.z;
			shift(n, &h);
			rotation_y(&h, steps[0].cosa, steps[0].sina);
			rotation_x(&h, steps[1].cosa, steps[1].sina);
			water[i].cp.p = sqrt(pow(h.x, 2) + pow(h.y, 2));
			water[i].cp.phi = atan(h.y / h.x);
			water[i].cp.z = h.z;
			if (h.z <= z_dir.z && h.z > 0) {
				water[i].for_calculate = 1;
			}
		}
	}
}

