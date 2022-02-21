# pandoc -s Readme.md -o Readme.pdf \
#     --title-prefix 'Diario de bordo' \
#     --normalize \
#     --smart \
#     --toc \
#     --latex-engine=`which xelatex`

pandoc -s Readme.md -o Readme.pdf \
    --title-prefix 'Diario de bordo' \
    --toc \
    --listings -H listings-setup.tex
