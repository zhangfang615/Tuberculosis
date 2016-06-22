__author__ = 'Fang Zhang'
__date__= '2016.5.3'
__email__= 'fza34@sfu.ca'
__function__ = 'check snp in reference'

import linecache
import sys

reference_path = sys.argv[1]
# reference_path = "/Users/Fang/Documents/Bioinformatics/test/polyTB_new/reference/H37Rv_AL123456.3.fasta"
SNP_vcf_path = sys.argv[2]
# SNP_vcf_path = "/Users/Fang/Documents/Bioinformatics/test/test_shell/results/ERR245646_resi.vcf"

SNP_file = file(SNP_vcf_path)

SNP_string=SNP_file.readline().strip()
while 1 and SNP_string:
    SNP_string=SNP_file.readline().strip()
    if SNP_string.startswith("Mycobacterium"):
        fields = SNP_string.split("\t")
        position=int(fields[1])
        nucleotide=fields[3]
        line = position/70
        line_position = position - line* 70
        if not line_position == 0:
            line += 2
            line_position -= 1
        else:
            line += 1
            line_position = 69
        theline = linecache.getline(reference_path, line)
        reference_nucleotide = theline[line_position]
        if nucleotide.upper() == reference_nucleotide.upper():
            print "SNP in "+ str(position) +" matches!"
        else:
            print "SNP "+nucleotide.upper()+" in "+ str(position) +" dismatches! "+reference_nucleotide.upper()+"\t"+SNP_string



