# Example of use
# sh -v prepare2commit.sh
# OR
# /usr/bin/time -v sh -v prepare2commit.sh

# tar -cvf - input/* | xz -v > input.tar.xz

tar -cvf - input/C174-PHASED.fa | xz -v > input/C174-PHASED.fa.xz
split -b 90m - input/C174-PHASED.fa.xz

tar -cvf - input/ResultadoFinal_C174-PHASED.gff | xz -v > input/ResultadoFinal_C174-PHASED.gff.fa.xz
tar -cvf - input/annot-genes/* | xz -v > input/annot-genes.fa.xz
tar -cvf - input/TEs/* | xz -v > input/TEs.fa.xz
