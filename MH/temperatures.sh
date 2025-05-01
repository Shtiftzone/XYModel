#!/bin/bash

TMIN=0.800
TMAX=1.100
N=61  # Iloœæ unikalnych temperatur
OUTFILE="/users/project1/pt01192/XYModels/MH/temperatures.txt"

rm -f "$OUTFILE"  # Usuwamy stary plik

for INDEX in $(seq 0 60); do
    T=$(awk -v min=$TMIN -v max=$TMAX -v idx=$INDEX -v total=$N 'BEGIN {print min + idx * (max - min) / (total - 1)}')
    TFORMATTED=$(printf "%.4f" "$T")  # Formatowanie temperatury
    printf "%s %d\n" "$TFORMATTED" "$INDEX" >> "$OUTFILE"
done

echo "Wygenerowano temperatures.txt z $(wc -l < $OUTFILE) temperaturami."