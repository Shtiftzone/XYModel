#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dirent.h>

#define L 64
#define NSPIN (L * L)
#define K 100
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int lut[] = {0, 1, 1, 0, 1, 0, 2, -1, 1, 2, 0, -1, 0, -1, -1, 0};

double mean_angle(double *theta) {
    double sumx = 0.0, sumy = 0.0;
    for (int i = 0; i < NSPIN; i++) {
        sumx += cos(theta[i]);
        sumy += sin(theta[i]);
    }
    return atan2(sumy, sumx);
}

int pbceuler(int *conf) {
    int euler = 0;
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

double average_euler(double *theta) {
    int conf[NSPIN];
    double mean = mean_angle(theta);
    double euler_sum = 0.0;
    for (int k = 0; k < K; k++) {
        double thresh = (double)k / K;
        for (int i = 0; i < NSPIN; i++) {
            double norm = fabs(atan2(sin(theta[i] - mean), cos(theta[i] - mean))) / M_PI;
            conf[i] = (norm >= thresh) ? 0 : 1;
        }
        int euler = pbceuler(conf);
        euler_sum += euler;
    }
    return euler_sum / K;
}

void process_file(const char *path, const char *filename) {
    char fullpath[512];
    snprintf(fullpath, sizeof(fullpath), "%s/%s", path, filename);
    FILE *fp = fopen(fullpath, "r");
    if (!fp) {
        perror("Error opening file");
        return;
    }

    int Nconf = 100;
    double theta[NSPIN];
    double sum = 0.0;
    double sum_sq = 0.0;

    for (int conf_idx = 0; conf_idx < Nconf; conf_idx++) {
        for (int i = 0; i < NSPIN; i++) {
            if (fscanf(fp, "%lf", &theta[i]) != 1) {
                fprintf(stderr, "Error reading file %s at configuration %d, index %d\n", filename, conf_idx, i);
                fclose(fp);
                return;
            }
        }
        double avg_euler = average_euler(theta);
        sum += avg_euler;
        sum_sq += avg_euler * avg_euler;
    }

    fclose(fp);

    double mean = sum / Nconf;
    double variance = (sum_sq / Nconf) - (mean * mean);

    printf("%s: Variance = %.8f\n", filename, variance);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s folder_with_txt_files\n", argv[0]);
        return 1;
    }

    const char *folder = argv[1];
    DIR *dir = opendir(folder);
    if (!dir) {
        perror("Error opening folder");
        return 1;
    }

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        if (strstr(entry->d_name, ".txt")) {
            process_file(folder, entry->d_name);
        }
    }

    closedir(dir);
    return 0;
}
