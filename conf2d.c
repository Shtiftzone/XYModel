//
//  conf2d.c
//
//  2-D cubic configurations for XY model with filtration-based topological descriptors
//
//  Functions:
//    mean_angle            - Computes average spin direction
//    pbceuler              - Euler number with periodic boundaries
//    write_me              - Magnetization and energy for XY model
//    write_filtered_topology - Average Euler, Betti, and face counts via filtration
//
#include <math.h>
#include <stdio.h>
#include <string.h>

#define NSPIN (L * L)
#define Z 4
#ifndef CELL
#define L2 (L * 2)
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int pbceuler(int *conf)
{
    int euler = 0;
    static int lut[] = {
#ifndef CELL
        0, 1, 1, 0, 1, 0, -2, -1, 1, -2, 0, -1, 0, -1, -1, 0
#else
        0, 1, 1, 0, 1, 0, 2, -1, 1, 2, 0, -1, 0, -1, -1, 0
#endif
    };
    for (int i = 0; i < L; i++) {
        int ip1 = (i + 1) % L;
        for (int j = 0; j < L; j++) {
            int jp1 = (j + 1) % L;
            euler += lut[conf[i + j * L] * 1 +
                         conf[ip1 + j * L] * 2 +
                         conf[i + jp1 * L] * 4 +
                         conf[ip1 + jp1 * L] * 8];
        }
    }
    return euler / 4;
}

double mean_angle(double *theta) {
    double sumx = 0.0, sumy = 0.0;
    for (int i = 0; i < NSPIN; i++) {
        sumx += cos(theta[i]);
        sumy += sin(theta[i]);
    }
    return atan2(sumy, sumx);
}

void write_me(double *theta, FILE *stream) {
    double mx = 0.0, my = 0.0, E = 0.0;
    for (int i = 0; i < L; i++) {
        int ip1 = (i + 1) % L;
        for (int j = 0; j < L; j++) {
            int jp1 = (j + 1) % L;
            int idx = i + j * L;
            int right = ip1 + j * L;
            int up = i + jp1 * L;

            mx += cos(theta[idx]);
            my += sin(theta[idx]);

            E -= cos(theta[idx] - theta[right]);
            E -= cos(theta[idx] - theta[up]);
        }
    }
    double M = sqrt(mx * mx + my * my);
    fprintf(stream, "%f,%f\n", M, E);
    fflush(stream);
}

void write_filtered_topology(double *theta, FILE *efile, FILE *bfile, FILE *ffile) {
    int conf[NSPIN];
    int data[NSPIN];
    int stack[NSPIN];
    double mean = mean_angle(theta);
    double euler_sum = 0.0, b0_sum = 0.0, b1_sum = 0.0, b2_sum = 0.0;
    double F0_sum = 0.0, F1_sum = 0.0, F2_sum = 0.0;
    /*double euler_values[K];
    double b0_values[K], b1_values[K], b2_values[K];
    double F0_values[K], F1_values[K], F2_values[K];*/

    for (int k = 0; k < K; k++) {
        double thresh = (double)k / K;
        for (int i = 0; i < NSPIN; i++) {
            double norm = fabs(atan2(sin(theta[i] - mean), cos(theta[i] - mean))) / M_PI;
            conf[i] = (norm >= thresh) ? 0 : 1 ;
        }

        int euler = pbceuler(conf);
        /*euler_values[k] = euler;*/
        euler_sum += euler;

        int b0 = 0, b2 = 1;
        memcpy(data, conf, NSPIN * sizeof(int));
        for (int i = 0; i < NSPIN; i++) {
            if (data[i] == 1) {
                b0++;
                int sp = 0;
                stack[sp++] = i;
                data[i] = 2;
                while (sp) {
                    int v = stack[--sp];
                    int x = v % L;
                    int y = v / L;
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            int nx = (x + dx + L) % L;
                            int ny = (y + dy + L) % L;
                            int nv = nx + ny * L;
                            if (data[nv] == 1) {
                                stack[sp++] = nv;
                                data[nv] = 2;
                            }
                        }
                    }
                }
            }
        }
        if (memchr(conf, 0, NSPIN * sizeof(int)) == NULL) b2 = 1; else b2 = 0;
        int b1 = b0 + b2 - euler;
        /*b0_values[k] = b0;
        b1_values[k] = b1;
        b2_values[k] = b2;*/
        b0_sum += b0; b1_sum += b1; b2_sum += b2;

        int F0 = 0, F1 = 0, F2 = 0;
        for (int i = 0; i < L; i++) {
            int ip1 = (i + 1) % L;
            for (int j = 0; j < L; j++) {
                int jp1 = (j + 1) % L;
                int quad = conf[i + j * L] * 1 +
                           conf[ip1 + j * L] * 2 +
                           conf[i + jp1 * L] * 4 +
                           conf[ip1 + jp1 * L] * 8;
                F0 += ((quad & 1) == 1);
                F1 += ((quad & 3) == 3) + ((quad & 5) == 5);
                F2 += (quad == 15);
            }
        }
        /*F0_values[k] = F0;
        F1_values[k] = F1;
        F2_values[k] = F2;*/
        F0_sum += F0; F1_sum += F1; F2_sum += F2;
    }

    /*for (int k = 0; k < K; k++) {
        fprintf(efile, "%d%c", (int)euler_values[k], (k < K - 1 ? ',' : '\n'));
        fprintf(bfile, "%d;%d;%d%c", (int)b0_values[k], (int)b1_values[k], (int)b2_values[k], (k < K - 1 ? ',' : '\n'));
        fprintf(ffile, "%d;%d;%d%c", (int)F0_values[k], (int)F1_values[k], (int)F2_values[k], (k < K - 1 ? ',' : '\n'));
    } */
    fprintf(efile, "%f\n", euler_sum / K);
    fprintf(bfile, "%f,%f,%f\n", b0_sum / K, b1_sum / K, b2_sum / K);
    fprintf(ffile, "%f,%f,%f\n", F0_sum / K, F1_sum / K, F2_sum / K);

    fflush(efile);
    fflush(bfile);
    fflush(ffile);
}