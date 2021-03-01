#include "protein_based_calculations.h"
/*
Main function is calculate_water_orientational_distribution
params:
    CA, C alpha coordinates vector
    CB, C beta coordinates vector
    N,  N coordinates vector
    HOHa, All hoh array
    HOHs, HOH coordinates array pointer
    distance: maximal distance from CA to O atom in A
return:
    int len : length of result coordinates array

*/

void rotation_x(vector* v, double cosa, double sina)
{
    double vy = (v->y * cosa) - (v->z * sina);
    double vz = (v->y * sina) + (v->z * cosa);
    v->y = check_if_null(vy);
    v->z = check_if_null(vz);
}

void rotation_y(vector* v, double cosa, double sina)
{
    /*x*/ double vx = (v->x * cosa) + (v->z * sina);
    /*z*/ double vz = ((-1) * v->x * sina) + (v->z * cosa);
    v->x = check_if_null(vx);
    v->z = check_if_null(vz);


}

void rotation_z(vector* v, double cosa, double sina)
{
    /*x*/ double vx = (v->x * cosa) - (v->y * sina);
    /*y*/ double vy = (v->x * sina) + (v->y * cosa);
    v->x = check_if_null(vx);
    v->y = check_if_null(vy);

}

void shift(vector* shift, vector* v)
{
    double x = v->x - shift->x;
    double y = v->y - shift->y;
    double z = v->z - shift->z;
    v->x = check_if_null(x);
    v->y = check_if_null(y);
    v->z = check_if_null(z);
}

double distance(vector* a, vector* b)
{
    return sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2) + pow(a->z - b->z, 2));
}

double check_if_null(double val)
{
    double check_val = 10 * pow(10, -10);
    double test_val = val;
    if (val < 0) {
        test_val = val * -1;
    }
    if (test_val < check_val) {
        return 0;
    }
    else {
        return val;
    }
}

double cos_a(vector* v1, vector* v2)
{
    double divident = (v1->x * v2->x) + (v1->y * v2->y) + (v1->z * v2->z);
    double divider = (sqrt(pow(v1->x, 2) + pow(v1->y, 2) + pow(v1->z, 2))) * sqrt((pow(v2->x, 2) + pow(v2->y, 2) + pow(v2->z, 2)));
    double result = divident / divider;
    return  check_if_null(result);
}

int calculate_water_orientational_distribution(vector* CA, vector* CB, vector* N, atom* HOHa, int HOHa_len, double dist, FILE* file_log)
{
    vector x_direction = { 1, 0, 0 }, y_direction = { 0, 1, 0 };

    //check which hohs in distance area
    int near_hoh_counter = 0;
    for (int i = 0; i < HOHa_len; i++) {
        if (distance(&HOHa[i], CA) < dist) {
            near_hoh_counter++;
        }
    }

    //Shift CA to coordinates origin. shift all other atoms to origin_shift

    shift(CA, CB);
    shift(CA, N);

    vector cb_oxy_proj = { CB->x, CB->y, 0 };
    step steps[3];
    ZeroMemory(steps, sizeof(step) * 3);

    double cosa, sina, alpha;
    double hemi_pi = M_PI / 2;
    // cosa CB - X
    cosa = cos_a(&cb_oxy_proj, &x_direction);
    if (CB->x < 0) {
        alpha = acos(cosa);
        cosa = cos((hemi_pi - alpha) + hemi_pi) * -1;
    }
    sina = sin_from_cos(cosa);
    if (cb_oxy_proj.y > 0) {
        sina = -sina;
    }
    steps[0].sina = sina;
    steps[0].cosa = cosa;
    rotation_z(CB, cosa, sina);
    rotation_z(N, cosa, sina);
    // cosa CB` - X
    cosa = cos_a(CB, &x_direction);
    if (CB->x < 0) {
        alpha = acos(cosa);
        cosa = cos((hemi_pi - alpha) + hemi_pi) * -1;
    }
    sina = sin_from_cos(cosa);
    if (CB->z < 0) {
        sina = -sina;
    }
    steps[1].sina = sina;
    steps[1].cosa = cosa;
    rotation_y(CB, cosa, sina);
    rotation_y(N, cosa, sina);
    if (N->z != 0) {
        // cosa N`` - X
        vector n_yoz_proj = { 0, N->y, N->z };
        cosa = cos_a(&n_yoz_proj, &y_direction);
        sina = sin_from_cos(cosa);
        if (N->z > 0) {
            sina = -sina;
        }

        steps[2].sina = sina;
        steps[2].cosa = cosa;
        rotation_x(N, cosa, sina);
    }
    else {
        steps[2].sina = 0;
        steps[2].cosa = 1;
    }

    if (CB->y != 0 || CB->z != 0) {
        print_vector(CB, "CB DATA");
        printf("CALCULATION ERROR:\tCB y or z not null\n");
        return -1;
    }

    if (N->z != 0) {
        print_vector(N, "N DATA");
        printf("CALCULATION ERROR:\tN->z not null\n");
        return -1;
    }

    for (int i = 0; i < HOHa_len; i++) {
        if (distance(&HOHa[i].c, CA) < dist) {

            atom* OH2, * H1, * H2;
            int lid;
            if (i == 0) {
                lid = 0;
            }
            else {
                lid = HOHa[i - 1].id;
            }
            int rid;
            if (i == HOHa_len - 1) {
                rid = HOHa[i + 1].id;
            }
            else {
                rid = 0;
            }
            if (strcmp(HOHa[i].name, "OH2") == 0 && &HOHa[i].for_calculate) {
                OH2 = &HOHa[i];
                H1 = &HOHa[i + 1];
                H2 = &HOHa[i + 2];
                i = i + 2;
            }
            else {
                continue;
            }
            if (
                check_if_null(OH2->c.x) == 0 ||
                check_if_null(OH2->c.y) == 0 ||
                check_if_null(OH2->c.z) == 0 ||
                check_if_null(H1->c.x) == 0 ||
                check_if_null(H1->c.y) == 0 ||
                check_if_null(H1->c.z) == 0 ||
                check_if_null(H2->c.x) == 0 ||
                check_if_null(H2->c.y) == 0 ||
                check_if_null(H2->c.z) == 0
                ) {
                continue;
            }
            shift(&CA, &OH2->c);
            shift(&CA, &H1->c);
            shift(&CA, &H2->c);

            rotation_z(&OH2->c, steps[0].cosa, steps[0].sina);
            rotation_y(&OH2->c, steps[1].cosa, steps[1].sina);
            rotation_x(&OH2->c, steps[2].cosa, steps[2].sina);

            rotation_z(&H1->c, steps[0].cosa, steps[0].sina);
            rotation_y(&H1->c, steps[1].cosa, steps[1].sina);
            rotation_x(&H1->c, steps[2].cosa, steps[2].sina);

            rotation_z(&H2->c, steps[0].cosa, steps[0].sina);
            rotation_y(&H2->c, steps[1].cosa, steps[1].sina);
            rotation_x(&H2->c, steps[2].cosa, steps[2].sina);

            shift(&OH2->c, &H1->c);
            shift(&OH2->c, &H2->c);

            double zenit1, zenit2, azimut1, azimut2, zo, ao, ro;

            vector zero = { 0, 0, 0 };
            ro = distance(&zero, &OH2->c);
            zo = acos(OH2->c.z / ro);
            ao = atan(OH2->c.y / OH2->c.x);

            zenit1 = acos(H1->c.z / (sqrt(pow(H1->c.x, 2) + pow(H1->c.y, 2) + pow(H1->c.z, 2))));
            zenit2 = acos(H2->c.z / (sqrt(pow(H2->c.x, 2) + pow(H2->c.y, 2) + pow(H2->c.z, 2))));

            //if (H1->c.x == 0) {
            //    azimut1 = 90;
            //}
            //if (H2->c.x == 0) {
            //    azimut2 = 90;
            //}

            azimut1 = atan(H1->c.y / H1->c.x);
            azimut2 = atan(H2->c.y / H2->c.x);

            //printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", zenit1, azimut1, zenit2, azimut2, ro, OH2->cp.phi);
            fprintf(file_log, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", zenit1, azimut1, zenit2, azimut2, ro, OH2->cp.phi);

        }
    }
    return 0;
}

double sin_from_cos(double cosa)
{
    double res = sqrt(1 - pow(cosa, 2));
    return check_if_null(res);
}

void print_vector(vector* a, const char* description)
{
    printf("VECTOR %s: x = %lf y = %lf z = %lf\n", description, a->x, a->y, a->z);
}
