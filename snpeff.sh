scp -r username@cluster:/path/nf-core_sarek/results/annotation/mutect2/sample_name/sample_name.mutect2.filtered_snpEff.ann.vcf.gz .

cat sample_name.mutect2.filtered_snpEff.ann.vcf \
    | ./scripts/vcfEffOnePerLine.pl \
    | java -jar SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].HGVS_P" "ANN[*].HGVS_C" "ANN[*].IMPACT" "ANN[*].BIOTYPE" " ( QUAL >= 30 )"  "(FILTER = 'PASS')" "ANN[*].GENE" "ANN[*].GENEID" > sample_name_NS_PASS.txt
