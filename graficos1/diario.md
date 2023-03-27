Genes, LncRNAs, miRNA, others ncRNAS (sem lncRNAs e sem miRNAs)

```bash
mkdir -p graf_CE
mkdir -p graf_CC
mkdir -p graf_CA
```

# Cromossomos

```bash
echo "chr;start;end" > CC_chr.csv
awk 'BEGIN {FS = ";";OFS = ";"} /^[0-9]+/ {print "chr"$1,$2,$3}' arquivos1/size_can.csv >> CC_chr.csv

echo "chr;start;end" > CE_chr.csv
awk 'BEGIN {FS = ";";OFS = ";"} /^[0-9]+/ {print "chr"$1,$2,$3}' arquivos1/size_eug.txt >> CE_chr.csv

echo "chr;start;end" > CA_C_chr.csv
awk 'BEGIN {FS = ";";OFS = ";"} /^sq00[0-9]{2}/ {print "chr"(substr($1,5)+0),$2,$3}' 4Alexandre/SG_C/size-sg_C.csv >> CA_C_chr.csv

echo "chr;start;end" > CA_E_chr.csv
awk 'BEGIN {FS = ";";OFS = ";"} /^sq00[0-9]{2}/ {print "chr"(substr($1,5)+0),$2,$3}' 4Alexandre/SG_E/size-sg_E.csv >> CA_E_chr.csv

wc -l arquivos1/size_eug.txt
wc -l arquivos1/size_can.csv
wc -l arquivos1/size-sg_?.csv
wc -l C*_chr.csv

mv -v CE_chr.csv graf_CE/CE_chr.csv
mv -v CC_chr.csv graf_CC/CC_chr.csv
mv -v CA_?_chr.csv graf_CA/
```

# Obter arquivos de Genes (1ª camada)

```bash
echo "chr;start;end;value1;color" > CC_genes.csv
awk 'BEGIN {FS = ";";OFS = ";"} /^[0-9]+/ {print "chr"$1,$2,$3,$4,"a"}' arquivos1/can_gene.csv >> CC_genes.csv

echo "chr;start;end;value1;color" > CE_genes.csv
awk 'BEGIN {FS = ";";OFS = ";"} /^[0-9]+/ {print "chr"$1,$2,$3,$4,"a"}' arquivos3-eug/CE1.5_v2_restitched_mueller.mRNA.gff.shiny.csv >> CE_genes.csv

echo "chr;start;end;value1;color" > CA_C_genes.csv
awk 'BEGIN {FS = ";";OFS = ";"} /^sq00[0-9]{2}/ {print "chr"(substr($1,5)+0),$2,$3,$4,"a"}' 4Alexandre/SG_C/sg_C_genes.csv >> CA_C_genes.csv

echo "chr;start;end;value1;color" > CA_E_genes.csv
# awk 'BEGIN {FS = ";";OFS = ";"} /^sq00[0-9]{2}/ {print "chr"(substr($1,5)+0),$2,$3,$4,"a"}' 4Alexandre/SG_E/sg_E_genes.csv >> CA_E_genes.csv
awk 'BEGIN {FS = ";";OFS = ";"} /^chr[0-9]{1,2}/ {print "chr"(substr($1,4)+0),$2,$3,$4,"a"}' arquivos3-SGE/CA_v0.6_v4_BTI_JS_edits_03.18.20.gff3.mRNA.shiny >> CA_E_genes.csv

wc -l arquivos1/*_gene.csv 4Alexandre/SG_?/sg_?_genes.csv
wc -l C*_genes.csv

mv -v CE_genes.csv graf_CE/CE_genes.csv
mv -v CC_genes.csv graf_CC/CC_genes.csv
mv -v CA_?_genes.csv graf_CA/
```

# Obter arquivos de TE (2ª camada)

Não é para fazer

# Obter arquivos de miRNA (4ª camada)

```bash
echo "chr;start;end;stack" > CC_miRNA.csv
awk 'BEGIN {OFS = ";"} /^Chr[0-9]{2}/ {if ($3=="miRNA") {print "chr"(substr($1,4)+0),$4,$5,"miRNA"}}' arquivos2/CC_ncRNAs_FINAL.gff3 >> CC_miRNA.csv

echo "chr;start;end;stack" > CE_miRNA.csv
awk 'BEGIN {OFS = ";"} /^chr(0[1-9]|[1-9][0-9])/ {if ($3=="miRNA") {print "chr"(substr($1,4)+0),$4,$5,"miRNA"}}' arquivos2/CE_ncRNAs_FINAL.gff3 >> CE_miRNA.csv

echo "chr;start;end;stack" > CA_C_miRNA.csv
awk 'BEGIN {OFS = ";"} /^sq00[0-9][13579]/ {if ($3=="miRNA") {print "chr"(substr($1,5)+0),$4,$5,"miRNA"}}' arquivos2/CA_ncRNAs_FINAL.gff3 >> CA_C_miRNA.csv

echo "chr;start;end;stack" > CA_E_miRNA.csv
awk 'BEGIN {OFS = ";"} /^sq00[0-9][02468]/ {if ($3=="miRNA") {print "chr"(substr($1,5)+0),$4,$5,"miRNA"}}' arquivos2/CA_ncRNAs_FINAL.gff3 >> CA_E_miRNA.csv
 
grep "^Chr" < arquivos2/CC_ncRNAs_FINAL.gff3 | grep "miRNA" | wc -l
grep "^chr" < arquivos2/CE_ncRNAs_FINAL.gff3 | grep "miRNA" | wc -l
grep "^sq00[0-9][13579]" < arquivos2/CA_ncRNAs_FINAL.gff3 | grep "miRNA" | wc -l
grep "^sq00[0-9][02468]" < arquivos2/CA_ncRNAs_FINAL.gff3 | grep "miRNA" | wc -l
wc -l C*_miRNA.csv

mv -v CC_miRNA.csv graf_CC/CC_miRNA.csv
mv -v CE_miRNA.csv graf_CE/CE_miRNA.csv
mv -v CA_?_miRNA.csv graf_CA/
```

# Obter arquivos de lncRNA (5ª camada)

```bash
echo "chr;start;end;stack" > CC_lncRNA.csv
awk 'BEGIN {OFS = ";"} /^Chr[0-9]{2}/ {print "chr"(substr($1,4)+0),$4,$5,"lncRNA"}' arquivos2/CC_LncRNAs_FINAL.gtf >> CC_lncRNA.csv

echo "chr;start;end;stack" > CE_lncRNA.csv
awk 'BEGIN {OFS = ";"} /^chr[0-9]{2}/ {print "chr"(substr($1,4)+0),$4,$5,"lncRNA"}' arquivos2/CE_LncRNAs_FINAL.gtf >> CE_lncRNA.csv

echo "chr;start;end;stack" > CA_C_lncRNA.csv
awk 'BEGIN {OFS = ";"} /^sq00[0-9][13579]/ {print "chr"(substr($1,5)+0),$4,$5,"lncRNA"}' arquivos2/CA_LncRNAs_FINAL.gtf >> CA_C_lncRNA.csv

echo "chr;start;end;stack" > CA_E_lncRNA.csv
awk 'BEGIN {OFS = ";"} /^sq00[0-9][02468]/ {print "chr"(substr($1,5)+0),$4,$5,"lncRNA"}' arquivos2/CA_LncRNAs_FINAL.gtf >> CA_E_lncRNA.csv

wc -l arquivos2/C?_LncRNAs_FINAL.gtf*
wc -l C*_lncRNA.csv

mv -v CC_lncRNA.csv graf_CC/CC_lncRNA.csv
mv -v CE_lncRNA.csv graf_CE/CE_lncRNA.csv
mv -v CA_?_lncRNA.csv graf_CA/
```

# Obter arquivos de ncRNA (3ª camada)

```bash
cp graf_CC/CC_lncRNA.csv CC_ncRNA.csv
cp graf_CE/CE_lncRNA.csv CE_ncRNA.csv
cp graf_CA/CA_C_lncRNA.csv CA_C_ncRNA.csv
cp graf_CA/CA_E_lncRNA.csv CA_E_ncRNA.csv
 
awk 'BEGIN {OFS = ";"} /^Chr[0-9]{2}/ {if ($3 ~ /RNA$/) {print "chr"(substr($1,4)+0),$4,$5,$3}}' arquivos2/CC_ncRNAs_FINAL.gff3 >> CC_ncRNA.csv
awk 'BEGIN {OFS = ";"} /^chr(0[1-9]|[1-9][0-9])/ {if ($3 ~ /RNA$/) {print "chr"(substr($1,4)+0),$4,$5,$3}}' arquivos2/CE_ncRNAs_FINAL.gff3 >> CE_ncRNA.csv
awk 'BEGIN {OFS = ";"} /^sq00[0-9][13579]/ {if ($3 ~ /RNA$/) {print "chr"(substr($1,5)+0),$4,$5,$3}}' arquivos2/CA_ncRNAs_FINAL.gff3 >> CA_C_ncRNA.csv
awk 'BEGIN {OFS = ";"} /^sq00[0-9][02468]/ {if ($3 ~ /RNA$/) {print "chr"(substr($1,5)+0),$4,$5,$3}}' arquivos2/CA_ncRNAs_FINAL.gff3 >> CA_E_ncRNA.csv

head -n 1 CC_ncRNA.csv > CC_ncRNA_v2.csv
awk 'BEGIN {FS = ";";OFS = ";"} {if (NR > 1) {print $1,$2,$3,"RNA"}}' CC_ncRNA.csv >> CC_ncRNA_v2.csv

head -n 1 CE_ncRNA.csv > CE_ncRNA_v2.csv
awk 'BEGIN {FS = ";";OFS = ";"} {if (NR > 1) {print $1,$2,$3,"RNA"}}' CE_ncRNA.csv >> CE_ncRNA_v2.csv

head -n 1 CA_C_ncRNA.csv > CA_C_ncRNA_v2.csv
awk 'BEGIN {FS = ";";OFS = ";"} {if (NR > 1) {print $1,$2,$3,"RNA"}}' CA_C_ncRNA.csv >> CA_C_ncRNA_v2.csv

head -n 1 CA_E_ncRNA.csv > CA_E_ncRNA_v2.csv
awk 'BEGIN {FS = ";";OFS = ";"} {if (NR > 1) {print $1,$2,$3,"RNA"}}' CA_E_ncRNA.csv >> CA_E_ncRNA_v2.csv

mv -v CC_ncRNA*.csv graf_CC/
mv -v CE_ncRNA*.csv graf_CE/
mv -v CA_?_ncRNA*.csv graf_CA/
```

# Obter arquivos de ncRNA - (lncRNA e miRNA)

```bash
echo "chr;start;end;stack" > CC_otherRNA.csv
awk 'BEGIN {OFS = ";"} /^Chr[0-9]{2}/ {if ($3 ~ /RNA$/ && $3 != "miRNA" && $3 != "lncRNA") {print "chr"(substr($1,4)+0),$4,$5,$3}}' arquivos2/CC_ncRNAs_FINAL.gff3 >> CC_otherRNA.csv

echo "chr;start;end;stack" > CE_otherRNA.csv
awk 'BEGIN {OFS = ";"} /^chr(0[1-9]|[1-9][0-9])/ {if ($3 ~ /RNA$/ && $3 != "miRNA" && $3 != "lncRNA") {print "chr"(substr($1,4)+0),$4,$5,$3}}' arquivos2/CE_ncRNAs_FINAL.gff3 >> CE_otherRNA.csv

echo "chr;start;end;stack" > CA_C_otherRNA.csv
awk 'BEGIN {OFS = ";"} /^sq00[0-9][13579]/ {if ($3 ~ /RNA$/ && $3 != "miRNA" && $3 != "lncRNA") {print "chr"(substr($1,5)+0),$4,$5,$3}}' arquivos2/CA_ncRNAs_FINAL.gff3 >> CA_C_otherRNA.csv

echo "chr;start;end;stack" > CA_E_otherRNA.csv
awk 'BEGIN {OFS = ";"} /^sq00[0-9][02468]/ {if ($3 ~ /RNA$/ && $3 != "miRNA" && $3 != "lncRNA") {print "chr"(substr($1,5)+0),$4,$5,$3}}' arquivos2/CA_ncRNAs_FINAL.gff3 >> CA_E_otherRNA.csv

head -n 1 CC_otherRNA.csv > CC_otherRNA_v2.csv
awk 'BEGIN {FS = ";";OFS = ";"} {if (NR > 1) {print $1,$2,$3,"RNA"}}' CC_otherRNA.csv >> CC_otherRNA_v2.csv

head -n 1 CE_otherRNA.csv > CE_otherRNA_v2.csv
awk 'BEGIN {FS = ";";OFS = ";"} {if (NR > 1) {print $1,$2,$3,"RNA"}}' CE_otherRNA.csv >> CE_otherRNA_v2.csv

head -n 1 CA_C_otherRNA.csv > CA_C_otherRNA_v2.csv
awk 'BEGIN {FS = ";";OFS = ";"} {if (NR > 1) {print $1,$2,$3,"RNA"}}' CA_C_otherRNA.csv >> CA_C_otherRNA_v2.csv

head -n 1 CA_E_otherRNA.csv > CA_E_otherRNA_v2.csv
awk 'BEGIN {FS = ";";OFS = ";"} {if (NR > 1) {print $1,$2,$3,"RNA"}}' CA_E_otherRNA.csv >> CA_E_otherRNA_v2.csv

mv -v CC_otherRNA*.csv graf_CC/
mv -v CE_otherRNA*.csv graf_CE/
mv -v CA_?_otherRNA*.csv graf_CA/
```
