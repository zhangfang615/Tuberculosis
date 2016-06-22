#!/bin/bash
# Shell script for counting SNPs on all samples
# Author: Fang Zhang
# Date: 2016.5.15
# E-mail: fza34@sfu.ca

# $var1 = path to list of samples
# $var2 = directory of results
# $var3 = path to vcf file of resistant SNPs

export PERL5LIB=/global/software/vcftools/vcftools_0.1.12b/perl:$PERL5LIB
export PATH=/global/software/htslib/htslib121/bin:$PATH
export PATH=/global/software/vcftools/vcftools_0.1.12b/bin:$PATH

# Reading parameters  ----------------------------------------
sample_list=${1}
results_path=${2}
count_file=${3}

cat ${sample_list} | while read sequence
do
    gatk_vcf=${results_path}"/"${sequence}"/"${sequence}"_gatk.vcf"
    gatk_count=`vcf-stats ${gatk_vcf} | grep -m1 'snp_count' | tr -s [:space:]|sed 's/.$//'|awk '{print $3}'`
    
    mpileup_vcf=${results_path}"/"${sequence}"/"${sequence}"_mpileup.vcf"
    mpileup_count=`vcf-stats ${mpileup_vcf} | grep -m1 'snp_count' | tr -s [:space:]|sed 's/.$//'|awk '{print $3}'`
    
    intersect_vcf=${results_path}"/"${sequence}"/"${sequence}"_intersect.vcf"
    intersect_count=`vcf-stats ${intersect_vcf} | grep -m1 'snp_count' | tr -s [:space:]|sed 's/.$//'|awk '{print $3}'`
    
#    freebayes_vcf=${results_path}"/"${sequence}"/"${sequence}"_freebayes.vcf"
#    freebayes_count=`vcf-stats ${freebayes_vcf} | grep -m1 'snp_count' | tr -s [:space:]|sed 's/.$//'|awk '{print $3}'`
    
    resi_vcf=${results_path}"/"${sequence}"/"${sequence}"_resi.vcf"
    resi_count=`vcf-stats ${resi_vcf} | grep -m1 'snp_count' | tr -s [:space:]|sed 's/.$//'|awk '{print $3}'`
    unresi_vcf=${results_path}"/"${sequence}"/"${sequence}"_noresi.vcf"
    unresi_count=`vcf-stats ${unresi_vcf} | grep -m1 'snp_count' | tr -s [:space:]|sed 's/.$//'|awk '{print $3}'`

    echo -e ${sequence}"\t"${gatk_count}"\t"${mpileup_count}"\t"${intersect_count}"\t"${resi_count}"\t"${unresi_count} >> ${count_file}
done

#  ./countSNPs.sh ../input/samples_list_0.txt ../test/test_results/results/ SNP_count_gatk.txt
