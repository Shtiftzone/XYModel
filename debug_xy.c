#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "conf2d.c"

#ifndef L
#define L 64
#endif

#ifndef K
#define K 100
#endif

#define NSPIN (L * L)

int main() {
    double theta[NSPIN];

    // Wczytaj theta z pliku
    FILE *f = fopen("/users/scratch1/shtiftzon/theta_input.txt", "r");
    if (!f) {
        perror("Nie mo¿na otworzyæ pliku theta_input.txt");
        return EXIT_FAILURE;
    }

    for (int i = 0; i < NSPIN; i++) {
        if (fscanf(f, "%lf", &theta[i]) != 1) {
            fprintf(stderr, "B³¹d wczytywania theta[%d]\n", i);
            return EXIT_FAILURE;
        }
    }
    fclose(f);

    // Wypisz wszystkie wartoœci charakterystyki Eulera dla filtracji
    int conf[NSPIN];
    double mean = mean_angle(theta);
    for (int k = 0; k < K; k++) {
        double thresh = (double)k / K;
        for (int i = 0; i < NSPIN; i++) {
            double diff = fabs(theta[i] - mean);
            if (diff > M_PI) diff = 2 * M_PI - diff;
            double norm = diff / M_PI;
            conf[i] = (norm >= thresh) ? 0 : 1;
        }

        int euler = pbceuler(conf);
        printf("%d%s", euler, (k < K - 1 ? "," : "\n"));
    }

    return 0;
}
