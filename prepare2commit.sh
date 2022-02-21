# Example of use
# sh -v prepare2commit.sh
# OR
# /usr/bin/time -v sh -v prepare2commit.sh

# tar -cvf - input/* | xz -v > input.tar.xz

tar -cvf - input/C174-PHASED.fa | xz -v > input/C174-PHASED.fa.xz
# split into smaller parts for github
split --verbose -b90M input/C174-PHASED.fa.xz input/C174-PHASED.fa.xz.
rm -v input/C174-PHASED.fa.xz
# to reassemble
# cat input/C174-PHASED.fa.xz.?? > input/C174-PHASED.fa.xz

tar -cvf - input/ResultadoFinal_C174-PHASED.gff | xz -v > input/ResultadoFinal_C174-PHASED.gff.fa.xz
tar -cvf - input/annot-genes/* | xz -v > input/annot-genes.fa.xz
tar -cvf - input/TEs/* | xz -v > input/TEs.fa.xz
