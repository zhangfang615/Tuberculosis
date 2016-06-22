#!/bin/bash
# Shell script for calling resistant SNPs on a single sample
# Author: Fang Zhang
# Date: 2016.5.13
# E-mail: fza34@sfu.ca

# $var1 = path to list of samples
# $var2 = directory of results
# $var3 = path to vcf file of resistant SNPs

# Reading parameters  ----------------------------------------
sample_list=${1}
results_path=${2}
ResiVCF=${3}
# position_file=${4}

cat ${sample_list} | while read sequence
do
    input_vcf=${results_path}"/"${sequence}"/"${sequence}"_intersect.vcf"
    echo ${input_vcf}
    input_gz=${input_vcf}".gz"
    echo ${input_gz}
    Sample_ResiVCF=${results_path}"/"${sequence}"/"${sequence}"_resi.vcf"
    bgzip -c ${input_vcf} > ${input_gz}
    tabix -p vcf ${input_gz}
    vcf-isec -o -n +2 ${input_gz} ${ResiVCF} > ${Sample_ResiVCF}
#    position=`grep Mycobacterium ${Sample_ResiVCF} | sed -n "2, 1p" | awk '{print $2}'`
#    echo ${position} >> ${position_file} 
done

# ./callResiSnps.sh ../input/samples_list_0.txt ../test/test_results/results/ ../input/RESI_SNPs/Resi-List-MasterV27_.vcf.gz
