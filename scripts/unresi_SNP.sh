#!/bin/bash
# Shell script for calling unresistant SNPs on a single sample
# Author: Fang Zhang
# Date: 2016.6.14
# E-mail: fza34@sfu.ca

# $var1 = path to list of samples
# $var2 = directory of results

export PERL5LIB=/global/software/vcftools/vcftools_0.1.12b/perl:$PERL5LIB
export PATH=/global/software/htslib/htslib121/bin:$PATH
export PATH=/global/software/vcftools/vcftools_0.1.12b/bin:$PATH

# Reading parameters  ----------------------------------------
sample_list=${1}
results_path=${2}

cat ${sample_list} | while read sequence
do
    intersect_vcf_gz=${results_path}"/"${sequence}"/"${sequence}"_intersect.vcf.gz"

    resi_vcf=${results_path}"/"${sequence}"/"${sequence}"_resi.vcf"
    
    resi_gz=${resi_vcf}".gz"
    
    unresi_vcf=${results_path}"/"${sequence}"/"${sequence}"_unresi.vcf"

    unresi_gz=${unresi_vcf}".gz"

    bgzip -c ${resi_vcf} > ${resi_gz}

    tabix -p vcf ${resi_gz}

    vcf-isec -c  ${intersect_vcf_gz} ${resi_gz} > ${unresi_vcf}

    bgzip -c ${unresi_vcf} > ${unresi_gz}

    tabix -p vcf ${unresi_gz}

    
done

