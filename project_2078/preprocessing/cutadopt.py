import os
import subprocess as sp
import glob
from multiprocessing import Pool

output_path = "/data/manke/group/rabbani/atac_project2078/trimmed_fastq"
def cutadapt(file):
	base_name = os.path.basename(file)
	dir_name = os.path.dirname(file)
	sample_name = base_name.split("_R1")[0]
	cmd = "module load cutadapt/2.10;"
	cmd += "cutadapt -j 10 -e 0.1 -q 16 -O 3 --trim-n --minimum-length 25 -A CTGTCTCTTATA -a CTGTCTCTTATA -o "+os.path.join(output_path, base_name)
	cmd += " -p "+os.path.join(output_path, sample_name+"_R2.fastq.gz")+" "
	cmd += file +" "+os.path.join(dir_name, sample_name+"_R2.fastq.gz")
	print(cmd)
	sp.check_output(cmd, shell = True)

with Pool(5) as p:
	p.map(cutadapt, glob.glob("/data/manke/group/rabbani/atac_project2078/fastq/*_R1.fastq.gz"))
