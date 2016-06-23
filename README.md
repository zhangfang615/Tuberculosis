Scripts used for Acquisition and transmission of TB drug resistance in China
========================================================================

documents
===============
    Experiments of determinism
    Experiments on sanity check

inputs
===============
	Beijing_lineage_list.txt: samples list for simulation
	shuffle_test.txt: shuffled samples list for determinism check
	RESI_SNPs: vcf file of resistant SNPs

reference
=============
	AL123456.2: H37Rv version AL123456.2
	AL123456.3: H37Rv version AL123456.3 (used in this research)
	TBSequence.fasta: H37Rv

results
============= 
    Beijing_lineage_SNP.xlsx: 
    SNP_count.txt		
    data_infomation_Beijing.xlsx
    determinism experiment.xlsx

simulation
=============
    TB_simulation.py: simulation

scripts
=============
    checkSNP.py: Python script for checking snp in reference
    countSNPs.sh: Shell script for counting SNPs on all samples
    callResiSnps.sh: Shell script for calling resistant SNPs on a single sample
    FileDownload.pbs: PBS script for file downloading
    FileDownload.py: Python script for file downloading
    jobSubmission.py: Python script for job submitting
    polyTB_Breezy.sh: Shell script for SNPs calling for one sample
    unresi_SNP.sh: Shell script for calling unresistant SNPs on a single sample