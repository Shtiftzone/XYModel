#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NRAND 512
#define Z 4
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Global variables
static double *theta;
static int L;
static int NSPIN;
static int irand = 0;
static double ran[NRAND];

// External random number functions
extern void ranlxd(double *r, int n);
extern void rlxd_init(int level, int seed);

// Initialize random seed
void sranlux(unsigned seed) {
    rlxd_init(1, seed);
}

// Random number
double dranlux(void) {
    if (irand == 0) {
        ranlxd(ran, NRAND);
        irand = NRAND - 1;
    }
    return ran[irand--];
}

// Initialize random configuration
void initialize_theta() {
    for (int i = 0; i < NSPIN; i++) {
        theta[i] = 2 * M_PI * dranlux();
    }
}

// Von Neumann neighbors (periodic)
int *neighbor(int i) {
    static int nn[Z];
    int p = 1;
    int q = 1 - L;
    int r = i;
    for (int j = 0; j < Z; j += 2) {
        nn[j] = (r + 1) % L ? i + p : i + q;
        nn[j + 1] = r % L ? i - p : i - q;
        p *= L;
        q *= L;
        r /= L;
    }
    return nn;
}

// Wolff update
void wolff_update(double T) {
    int seed = (int)(dranlux() * NSPIN);
    
    double angle_r = 2.0 * M_PI * dranlux(); // random reflection direction
    double rx = cos(angle_r);
    double ry = sin(angle_r);

    int visited[NSPIN];
    memset(visited, 0, NSPIN * sizeof(int)); // <- poprawka

    int stack[NSPIN];
    int sp = 0;
    stack[sp++] = seed;
    visited[seed] = 1;

    while (sp) {
        int current = stack[--sp];
        int *nn = neighbor(current);

        double sx = cos(theta[current]);
        double sy = sin(theta[current]);

        for (int j = 0; j < Z; j++) {
            int n = nn[j];
            if (visited[n]) continue;

            double snx = cos(theta[n]);
            double sny = sin(theta[n]);
            double delta_E = -2.0 * (rx * sx + ry * sy) * (rx * snx + ry * sny);
            double p_add = 1.0 - exp(delta_E / T);

            if (dranlux() < p_add) {
                visited[n] = 1;
                stack[sp++] = n;
            }
        }

        // Reflect spin across direction r
        double proj = rx * cos(theta[current]) + ry * sin(theta[current]);
        double refl_x = cos(theta[current]) - 2.0 * proj * rx;
        double refl_y = sin(theta[current]) - 2.0 * proj * ry;
        theta[current] = atan2(refl_y, refl_x);
    }
}

// Measure magnetization and energy
void measure(double *M, double *E) {
    double mx = 0.0, my = 0.0, energy = 0.0;
    for (int i = 0; i < L; i++) {
        int ip1 = (i + 1) % L;
        for (int j = 0; j < L; j++) {
            int jp1 = (j + 1) % L;
            int idx = i + j * L;
            int right = ip1 + j * L;
            int up = i + jp1 * L;

            mx += cos(theta[idx]);
            my += sin(theta[idx]);
            energy -= cos(theta[idx] - theta[right]);
            energy -= cos(theta[idx] - theta[up]);
        }
    }
    *M = sqrt(mx * mx + my * my) / NSPIN;
    *E = energy / NSPIN;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s L steps_per_conf\n", argv[0]);
        return 1;
    }

    L = atoi(argv[1]);
    NSPIN = L * L;
    int steps_per_conf = atoi(argv[2]);

    int Nconf = 10;
    unsigned seed = 123456;
    sranlux(seed);

    theta = malloc(NSPIN * sizeof *theta);
    if (!theta) {
        fprintf(stderr, "Couldn't allocate memory.\n");
        return 1;
    }

    double T_start = 0.88;
    double T_end = 1.1;
    double T_step = 0.01;

    for (double T = T_start; T <= T_end + 1e-6; T += T_step) {
        char filename[512];
        snprintf(filename, sizeof(filename), "./thermalization_L%d_T%.3f.txt", L, T);
        FILE *out = fopen(filename, "w");
        if (!out) {
            perror("Couldn't open output file");
            printf("Attempted to open: %s\n", filename);
            free(theta);
            return 1;
        }

        fprintf(out, "# Step Magnetization Energy\n");

        for (int conf = 0; conf < Nconf; conf++) {
            initialize_theta();
            fprintf(out, "# Configuration %d\n", conf);
            for (int step = 0; step < steps_per_conf; step++) {
                wolff_update(T);
                double M, E;
                measure(&M, &E);
                fprintf(out, "%d %.8f %.8f\n", step, M, E);
            }
        }

        fclose(out);
        printf("Finished T=%.3f and wrote to %s\n", T, filename);
    }

    free(theta);
    return 0;
}
