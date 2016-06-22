import subprocess
import sys
import time

File_path = "/home/zhf615/TB_test/MALAWI/input/Beijing_lineage_list_test.txt"
input_file = open(File_path)
sequence = '1'
while 1 and sequence:
	sequence = input_file.readline().strip()
	p = subprocess.Popen(["qsub", "-d", "/home/zhf615/TB_test/MALAWI/test",  "-l", "mem=160000mb", "-N", sequence, "-v", "var1=shuffled,var2=/home/zhf615/TB_test/MALAWI/reference/AL123456.3,var3=TBSequence,var4="+sequence+",var5=/global/scratch/zhf615/MALAWI/data/fastq/,var6=/home/zhf615/TB_test/MALAWI/input/RESI_SNPs/Resi-List-MasterV27_.vcf.gz", "./polyTB_Breezy.sh"], stdout=subprocess.PIPE)
	output, err = p.communicate()
	print "*** Submit job ***\n", output
#	time.sleep(20)
#        p = subprocess.Popen(["qsub", "-d", "/home/zhf615/TB_test/MALAWI/test", "-l", "mem=3000mb", "-N", sequence, "-v", "var1=normal,var2=/home/zhf615/TB_test/MALAWI/reference/AL123456.3,var3=TBSequence,var4="+sequence+",var5=/global/scratch/zhf615/MALAWI/data/fastq/,var6=/home/zhf615/TB_test/MALAWI/input/RESI_SNPs/Resi-List-MasterV27_.vcf.gz", "./polyTB_Breezy.sh"], stdout=subprocess.PIPE)
#        output, err = p.communicate()
#        print "*** Submit job ***\n", output
        time.sleep(50)
input_file.close()
# qsub -d /home/zhf615/TB_test/ -N polyTB -v var1='sample',var2='/home/zhf615/TB_test',var3='TBSequence',var4='ERR124635',var5='./samples' ./polyTB_Breezy.sh


