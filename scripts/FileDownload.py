__author__ = 'Fang'

import sys
import urllib

doc_path = sys.argv[1]
result_path = sys.argv[2]
# doc_path ="/Users/Fang/Documents/Bioinformatics/test/test_shell/filelist.txt"
# result_path="/Users/Fang/Documents/Bioinformatics/test/test_shell/result.txt"
infile = file(doc_path)
output_file = file(result_path, "w+")
while 1:
    input_line = infile.readline().strip()
    if not input_line:
        break
    surfix = input_line[0:6]
    url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+surfix+'/' + input_line +'/'+input_line+'_1.fastq.gz'
    # print "downloading : " + input_line
    try:
        urllib.urlretrieve(url, input_line+"_1.fastq.gz")
    except:
        print input_line+"_1.fastq.gz downloaded failed."
        output_file.writelines(input_line+"_1.fastq.gz downloaded failed.\r\n")
    url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/'+surfix+'/' + input_line +'/'+input_line+'_2.fastq.gz'
    # print "downloading : " + input_line
    try:
        urllib.urlretrieve(url, input_line+"_2.fastq.gz")
    except:
        print input_line+"_2.fastq.gz downloaded failed."
        output_file.writelines(input_line+"_2.fastq.gz downloaded failed.\r\n")
infile.close()
output_file.close()
