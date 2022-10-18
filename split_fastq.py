import argparse
import os 
import sys
import gzip
import re

## extract flowcell ID and lane as read group
def parse_header(x):
    temp=x.split(":")
    RG=".".join([temp[2],"L"+temp[3]])
    return(RG)

def get_readname(x):
    temp=x.split(" ")
    return(temp[0])

def get_readN(x):
    return(x.split(" ")[1][0])

def read_Nline(f,N=3):
    return("".join([f.readline() for x in range(N)]))

def writeFastq(reads,out,prefix,r):
    for i in reads.keys():
        fileout=f"{out}/{prefix}_{i}_R{r}_001.fastq.gz"
        out_write=gzip.open(fileout,"wt")
        out_write.writelines("".join(reads[i]))
        out_write.close()

def readFastq(fqFile,prefix):
    counter=0
    outfile={}
    while fqFile[0]: ## using the first file to go through all lines 
        if((counter % 1000000)==0):
                print(f"Processed {counter*len(fqFile)} lines...",flush=True)
        lines=[x.readline() for x in fqFile]
        if(lines[0]==""): # EOF
            break
        check_head=[x[0]=="@" for x in lines]
        if(sum(check_head)==len(lines)): ## this is a header line
            RG=[parse_header(x) for x in lines]
            read_names=[get_readname(x) for x in lines]
            read_N=[get_readN(x) for x in lines]
            if(len(set(read_names))!=1): ## the read groups are different between files
                    sys.stderr.write(f"Read names do not match at line {counter}.\n")
                    exit(3)
            ### read the next three lines
            seqs=[read_Nline(x,3) for x in fqFile] #read the next three lines
            out=["".join(x) for x in zip(lines,seqs)] ## concatenating header and the next three sequences
        elif(sum(check_head)!=0): ## headers are not matched in two files, throw an error 
            sys.stderr.write(f"Line orders don't seem to match in two FASTQ files at line {counter}.\n")
            exit(2)
        for rg,data,R in zip(RG,out,read_N): ## more efficient to store a list of string and joint than cat
            fileout=f"{prefix}_{rg}_R{R}_001.fastq.gz"
            if(fileout in outfile.keys()):
                outfile[fileout].write(data)
            else:
                if(os.path.exists(fileout)): ## if file exists, overwrites it
                    outfile[fileout]=gzip.open(fileout,"wt")
                else:
                    outfile[fileout]=gzip.open(fileout,"at")
                outfile[fileout].write(data)
            [x.close for x in outfile.values()]
        counter+=4


def main():
        parser= argparse.ArgumentParser(description="Split a FASTQ file into multiple files based on flowcell ID and lane.\nWill create a new directory with all new FASTQs within it.",formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("fastq",nargs="+",help="A FASTQ file. If paired-end data, the R2 file can also be provided.")
        parser.add_argument("--out",dest="out",help="Specify the output directory. Writes to a directory named temp by default. Note some existing FASTQs in the --out directory might be overwritten.")
        parser.add_argument("--prefix",dest="prefix",help="Specify the prefix of for new FASTQ files. Default is to use the prefix of provided FASTQs")
        args = parser.parse_args()
        if(len(args.fastq)>2):
            sys.stderr.write("More than 2 FASTQ files were provided. Please only provide one or two files for single-end or paired-end reads, respectively\n")
            exit(1)

        filename=[os.path.basename(x) for x in args.fastq]  ## get fastq file name

        if(args.prefix==None): ## using the first file prefix as prefix
            args.prefix=re.subn("_R1_001.fastq.gz|_R2_001.fastq.gz","",filename[0])[0]
        if(args.out==None):
            args.out="temp"
        os.system(f'mkdir -p {args.out}')  ## create output directory
        
        fqFile=[gzip.open(x,"rt") for x in args.fastq]
        reads=readFastq(fqFile,f"{args.out}/{args.prefix}") ## reads fastq
        [x.close() for x in fqFile]
        

if __name__ == "__main__":
    main()
