
import multiprocessing
from functools import partial
from multiprocessing import Process, Value,  freeze_support,  current_process,get_context
import sys
import subprocess
import os
import sys
import pysam
import pybedtools
import pysam.bcftools as bcftools
import numpy as np
import pandas as pd
from cyvcf2 import VCF
import argparse
import scripts.cnvest as cnv
import scripts.cnvtab as ct
import re
import glob
import pysam
import time
from itertools import repeat
from multiprocessing import set_start_method
def prep_fun(samples,bamfolder,bedfile,workdir,threads):
    #job_pool
    freeze_support()
    pool = multiprocessing.Pool(int(threads))
    
    global didx
    didx=0
    argms=[bamfolder,bedfile,workdir,threads]
    print(argms)
    df_cov=pool.starmap(depth_cal,zip(repeat(argms),samples))
    indx=[]
    for i in range(len(samples)):
        dfs_cov.append(df_cov[i])#globals()[f'df_cov{i}'])
        indx.append(i)  
    df_sam,dp_prop=cnv.cnv_cal(samples,dfs_cov,indx)
    for samidx in range(len(samples)):
        print(samples[samidx],df_sam.iloc[samidx])
        outfilename=workdir+"/output/report/"+str(samples[samidx])+"_depth_cal.txt"
        df_sam.iloc[samidx].to_csv(outfilename, index=False,sep="\t")
        ct.make_tab(df_sam.iloc[samidx],samples[samidx],workdir)


def depth_cal(ps,sample):
    print(sample)
    (bamfolder,bedfile,workdir,threads) =ps

    global dfs_cov
    for bamf in os.listdir(bamfolder):
        if bamf.endswith('.bam'):
            bam=bamfolder+"/%s-ready.bam" %sample
        if bamf.endswith('.cram'):
            bam=bamfolder+"/%s.cram" %sample
    if os.path.isfile(bam) == False:
        raise Exception("Bam file does not exist")
    process = current_process()
    # report the name of the process
    pn=process.name
    didx=0
    print("processno",get_context('fork'))
    pidx=pn.replace("ForkPoolWorker-", "")
    pidx= didx+ int(pidx) -1

    print("processing " + bam)
    bed_df= pd.read_csv(bedfile,sep='\t',usecols=[0,1,2,3,4], names=["chr","start","end","gene","strand"])
    #Directory and file list
    parentdir=workdir+'/output/'
    reportdir=workdir+'/output/report'
    outdir=workdir+'/output/vcf'
    tmpdir=workdir+'/output/tmp'
    if not os.path.exists(parentdir):
        os.mkdir(parentdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(reportdir):
        os.mkdir(reportdir)
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    
    tmpfile=tmpdir+'/'+str(sample)+'.'+'depthbcf_tmp.txt'
    tmpfile1=tmpdir+'/'+str(sample)+'.'+'depth_tmp.txt'
    tmpfile2=tmpdir+'/'+str(sample)+'.'+'tmp.bcf'
    outfile1=outdir+'/'+str(sample)+'.'+'depth.vcf'
    outfile=reportdir+'/'+str(sample)+'.'+'out.txt'
    
    if os.path.exists(outfile1):
        os.remove(outfile1)
        print("Vcf file exists, overwritten")
    if os.path.isfile(bam) == False:
        raise Exception("Bam file does not exist!")
    if os.path.isfile(bedfile) == False:
        raise Exception("Bed file does not exist!")
    
    #Vcf 
    df_name = pd.DataFrame(columns=['chr','Pos','Gene','DP', 'GT','MQ','DP4','start','end'])
    #pysam.depth("-a", bam,"-o",outfile,"-r",bedfile,catch_stdout=False)
    for ib in range(len(bed_df)) :
        bchr=bed_df.loc[ib, "chr"]
        bstart=bed_df.loc[ib, "start"]
        bend=bed_df.loc[ib, "end"]
        bgene=bed_df.loc[ib, "gene"]
        bstrand=bed_df.loc[ib, "strand"]
        bpos=bchr+':'+str(bstart)+'-'+str(bend)
        print("bpos",bpos)
        
    #,"--no-BAQ", "-d","100000000","-BQ","0","-Q", "0",
    bcftools.mpileup("--count-orphans","-f","data/GRCh38_full_analysis_set_plus_decoy_hla.fa", bam,"-R",bedfile, "-o",tmpfile2,catch_stdout=False)
    bcftools.call("-c" , "-o",tmpfile, tmpfile2,catch_stdout=False) #,"--insert-missed"
    with open(tmpfile,'r') as fp: 
        data=fp.read()
    with open(outfile1,'a') as op:
        op.write(data)
            
        
    a = pybedtools.BedTool.from_dataframe(bed_df)
    b = pybedtools.BedTool(outfile1)
    c= a.intersect(b,wb=True)
    c_df=pybedtools.BedTool.to_dataframe(c,disable_auto_names=True, header=None)
    #print("c_df",c_df,sample)
    df_name['chr'] = c_df.iloc[:,5]
    df_name['Pos'] = c_df.iloc[:,6]
    df_name['DP']= c_df.iloc[:,12].str.extract(r'DP=([0-9]+);', expand = True)
    df_name['GT']=c_df.iloc[:,14].str.extract(r'(\d{0,2}\/\d{0,2})\:\d{0,2}', expand = True)
    df_name['MQ']=c_df.iloc[:,12].str.extract(r'MQ=([0-9]+);', expand = True)
    df_name['Gene'] = c_df.iloc[:,3]
    df_name['DP4']= c_df.iloc[:,12].str.extract(r'DP4=([0-9,]+);', expand = True)
    df_name['start'] = c_df.iloc[:,1]
    df_name['end'] = c_df.iloc[:,2]
    
    df_name.drop_duplicates(subset=['chr', 'Pos'],keep='first',inplace=True) #To remove duplicate values if any
    df_name=df_name.reset_index(drop=True)

    ##Remove temporary files
    if os.path.exists(tmpfile):
        os.remove(tmpfile)
    if os.path.exists(tmpfile1):
        os.remove(tmpfile1)
    if os.path.exists(tmpfile2):
        os.remove(tmpfile2)
    #print(df_name)
    ##Filter low coverage positions
    df_name['DP'] = pd.to_numeric(df_name['DP'])
    #df_name.drop(df_name[df_name['DP'] < 20].index, inplace = True)
    #print(sample,df_name)

    #print(dfs_cov)

    return(df_name)

def main():
    indx = []
    global dfs_cov 
    dfs_cov=[] #coverage dataframe
    
    didx=0
    #set_start_method('spawn')
    

    ##Input arguments
    parser = argparse.ArgumentParser(description='CNVexpo ')
    parser.add_argument('--infile', required=True,help='sample list')
    parser.add_argument('--bedfile', required=True,help='Target bed file')
    parser.add_argument('--bamfolder', required=True,help='location of the bam files.')
    parser.add_argument('--workdir', required=True,help='location of the work directory.')
    parser.add_argument('--threads', required=True,help='number of cpus.')
    args = parser.parse_args()
    fname=args.workdir+'/input/'+args.infile
    bname=args.workdir+'/input/'+args.bedfile
    bamfolder=args.bamfolder
    workdir=args.workdir

    #Multiprocessing
    threads= args.threads
    lines = [line.strip() for line in open(fname, 'r')]

    prep_fun(lines,bamfolder,bname,workdir,threads )


if __name__ == "__main__":
    main()