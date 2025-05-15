//
//  xy.c
//
//  XY model with Wolff updates and filtration-based topology export
//
//  Parameters:
//    T                 - Temperature
//    Nconf             - Number of configurations
//
//  Inputs:
//    Sfile             - Seed
//    Ndecorr           - Number of decorrelation steps (optional)
//
//  Outputs:
//    mefile            - Magnetization and energy
//    efile             - Average Euler number from filtration
//    bfile             - Average Betti numbers from filtration
//    ffile             - Average face counts from filtration
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ranlux.h"

#define NRAND 360
#define K 100

#ifndef DIM
#include "conftr.c"
#elif DIM == 2
#include "conf2d.c"
#elif DIM == 3
#include "conf3d.c"
#elif DIM == 4
#include "conf4d.c"
#else
#error Unsupported configuration.
#endif

static double *theta;

void sranlux(const char *filename)
{
    unsigned seed;
    FILE *stream = fopen(filename, "r");
    if (stream && fscanf(stream, "%u\n", &seed) == 1) {
        fclose(stream);
        rlxd_init(1, seed);
    } else {
        stream = fopen("/dev/urandom", "r");
        if (stream == NULL || fread(&seed, sizeof seed, 1, stream) != 1) {
            fprintf(stderr, "Couldn't read from /dev/urandom\n");
            exit(EXIT_FAILURE);
        }
        fclose(stream);
        seed &= 0x3fffffff;
        rlxd_init(1, seed);
        stream = fopen(filename, "w");
        if (stream == NULL || fprintf(stream, "%u\n", seed) < 0) {
            fprintf(stderr, "Couldn't write to %s\n", filename);
            exit(EXIT_FAILURE);
        }
        fclose(stream);
    }
}

double dranlux(void)
{
    static int irand = 0;
    static double ran[NRAND];
    if (irand == 0) {
        ranlxd(ran, NRAND);
        irand = NRAND - 1;
    }
    return ran[irand--];
}

FILE *file_notna(const char *filename)
{
    return strcmp(filename, "NA") ? fopen(filename, "w") : NULL;
}

#ifdef DIM
int *neighbor(int i)
{
    // Von Neumann neighborhood on a cubic lattice
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
#endif

void wolff_update(double T) {
    // Select a random spin as the seed
    int seed = (int)(dranlux() * NSPIN);

    // Choose a random reflection direction (unit vector)
    double refl_angle = 2 * M_PI * dranlux();
    double rx = cos(refl_angle);
    double ry = sin(refl_angle);

    int visited[NSPIN] = {0};
    int stack[NSPIN];
    int sp = 0;
    stack[sp++] = seed;
    visited[seed] = 1;

    while (sp) {
        int current = stack[--sp];
        int *nn = neighbor(current);

        for (int j = 0; j < Z; j++) {
            int n = nn[j];
            if (visited[n]) continue;

            // Project spins onto the reflection direction
            double proj = rx * cos(theta[n]) + ry * sin(theta[n]);
            double proj_current = rx * cos(theta[current]) + ry * sin(theta[current]);

            // Calculate the bond probability
            double delta = proj * proj_current;
            if (delta > 0) { // Only if the projection is positive
                double p_add = 1.0 - exp(-2.0 * delta / T);
                if (dranlux() < p_add) {
                    visited[n] = 1;
                    stack[sp++] = n;
                }
            }
        }

        // Reflect the current spin across the reflection direction
        double sx = cos(theta[current]);
        double sy = sin(theta[current]);
        double dot = sx * rx + sy * ry;
        double refl_x = sx - 2.0 * dot * rx;
        double refl_y = sy - 2.0 * dot * ry;
        theta[current] = atan2(refl_y, refl_x);
    }
}

int main(int argc, char *argv[]) {
    int Ntherm = (int)(0.5 * L * L);
    int Ndecorr = 50;
    if (argc < 8) {
        fprintf(stderr,
                "%s T Nconf Sfile mefile bfile efile ffile [Ndecorr]\n",
                argv[0]);
        return EXIT_FAILURE;
    }
    theta = malloc(NSPIN * sizeof *theta);
    if (!theta) {
        fprintf(stderr, "Couldn't allocate memory\n");
        return EXIT_FAILURE;
    }
    double T = strtod(argv[1], NULL);
    int Nconf = (int) strtol(argv[2], NULL, 10);
    sranlux(argv[3]);
    FILE *Mefile = file_notna(argv[4]);
    FILE *Bfile = file_notna(argv[5]);
    FILE *Efile = file_notna(argv[6]);
    FILE *Ffile = file_notna(argv[7]);
    if (argc == 9)
        Ndecorr = (int) strtol(argv[8], NULL, 10);

    for (int i = 0; i < NSPIN; i++)
        theta[i] = 2 * M_PI * dranlux();
    for (int i = 0; i < Ntherm; i++)
        wolff_update(T);

    while (Nconf--) {
        if (Mefile)
            write_me(theta, Mefile);
        if (Efile && Bfile && Ffile)
            write_filtered_topology(theta, Efile, Bfile, Ffile);

        for (int i = 0; i < Ndecorr; i++)
            wolff_update(T);
    }
    return 0;
}
