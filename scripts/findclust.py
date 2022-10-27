import multiprocessing
from functools import partial
from multiprocessing import Process, Value,  freeze_support,  current_process,get_context
from itertools import repeat
from multiprocessing import set_start_method
from sklearn.cluster import MeanShift
from sklearn.cluster import  estimate_bandwidth
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut
from functools import reduce
from sklearn import preprocessing
import numpy as np
import os
import re
import pysam
import pybedtools
import pysam.bcftools as bcftools
import numpy as np
import pandas as pd
import argparse
from cyvcf2 import VCF
pd.set_option('use_inf_as_na',True)

def prep_fun(samples,bamfolder,workdir,threads):
    freeze_support()
    pool = multiprocessing.Pool(int(threads))
    #job_pool
    global didx
    didx=0
    argms=[bamfolder,workdir,threads]
    df_cov=pool.starmap(depth_cal,zip(repeat(argms),samples))
    tmpdir = workdir+'/output/tmp'
    outname = 'pd_tmp.txt'
    fullname = os.path.join(tmpdir, outname) 
    name=[]
    indx=[]
    dp_df=pd.DataFrame()
    dp_dfn=pd.DataFrame()
    for i in range(len(samples)):
        indx.append(i)
        if i==0:
            dp_df['chr']=df_cov[i]['chr']
            dp_df['Pos']=df_cov[i]['Pos']
            dp_df['Gene']=df_cov[i]['Gene']
            dp_df['chr'].astype('str')
            dp_df['Pos'].astype('int')
            dp_df['Gene'].astype('str')
            dp_df.to_csv(fullname,sep='\t')
        dp_df.drop_duplicates(keep=False, inplace=True)
        dp_dfn=pd.DataFrame({'chr':df_cov[i]['chr'],'Pos':df_cov[i]['Pos'],samples[i]:df_cov[i]['DP']})
        dp_dfn.iloc[:,3:] = dp_dfn.iloc[:,3:].astype(int)
        dp_df=pd.read_csv(fullname,sep='\t', index_col=None)
        dp_df.drop_duplicates(keep=False)
        dp_df.drop(dp_df.columns[dp_df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
        dp_df=dp_df.fillna(0)
        dp_df['chr'].astype('str')
        dp_df['Pos'].astype('int')
        dp_df=pd.merge(dp_df,dp_dfn,on=['chr','Pos'])
        dp_df.iloc[:,3:] = dp_df.iloc[:,3:].astype(int)
        dp_df.rename(columns={ dp_df.columns[i+3]: samples[i] }, inplace = True)
        dp_df[samples[i]].astype('int')
        dp_df.drop(dp_df.columns[dp_df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
        dp_df.drop_duplicates(keep=False, inplace=True)
        dp_df.info(memory_usage="deep")
        f = open(fullname, 'r+')
        f.truncate(0)  
        dp_df.to_csv(fullname,sep='\t', index=False)
    headers=pd.DataFrame()
    headers =  dp_df.iloc[:,3:].groupby(dp_df.iloc[:,2]).sum()
    headers = headers.rename_axis('Gene').reset_index()
    dp_df.replace([np.inf, -np.inf], np.nan, inplace=True) 
    dp_df = dp_df.dropna()
    dp_df = dp_df.reset_index( drop=True    )
    print("dp_df",dp_df.iloc[:,3:].astype(float))
    dp_df.iloc[:,3:]= dp_df.iloc[:,3:].astype(float)
    dp_df.dropna(axis=0)  
    sam_sel(dp_df,workdir)

def depth_cal(ps, sample):
    (bamfolder,workdir,threads)=ps
    for bamf in os.listdir(bamfolder):
        if bamf.endswith('.bam'):
            bam=bamfolder+"/%s-ready.bam" %sample
        if bamf.endswith('.cram'):
            bam=bamfolder+"/%s.cram" %sample
    if os.path.isfile(bam) == False:
        raise Exception("Bam file does not exist")

    bedfile="data/control.bed"
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

    for ib in range(len(bed_df)) :
        bchr=bed_df.loc[ib, "chr"]
        bstart=bed_df.loc[ib, "start"]
        bend=bed_df.loc[ib, "end"]
        bgene=bed_df.loc[ib, "gene"]
        bstrand=bed_df.loc[ib, "strand"]
        bpos=bchr+':'+str(bstart)+'-'+str(bend)
        print("bpos",bpos)


    bcftools.mpileup("--count-orphans","-f","data/GRCh38_full_analysis_set_plus_decoy_hla.fa",bam,"-R",bedfile, "-o",tmpfile2,catch_stdout=False)
    bcftools.call("-c" , "-o",tmpfile, tmpfile2,catch_stdout=False) 
    with open(tmpfile,'r') as fp: 
        data=fp.read()
    with open(outfile1,'a') as op:
        op.write(data)
    df_name=pd.DataFrame()
    a = pybedtools.BedTool.from_dataframe(bed_df)
    b = pybedtools.BedTool(outfile1)
    c= a.intersect(b,wb=True)
    c_df=pybedtools.BedTool.to_dataframe(c,disable_auto_names=True, header=None)
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
    df_name['DP'] = pd.to_numeric(df_name['DP'])
    return(df_name)
    
def sam_sel(dp_df,workdir):
    lt=len(dp_df.columns)
    print("sample",dp_df)
    dp_df.iloc[:,3:] = preprocessing.scale(dp_df.iloc[:,3:])
    newdp=dp_df.iloc[:,3:].mean()
    newdp=np.array(newdp).astype(float)
    bw = estimate_bandwidth(newdp[:,None], quantile=0.3)
    print(newdp)
    kernel='gaussian'
    ms = MeanShift(bandwidth=bw, bin_seeding=True)
    ms.fit(newdp[:,None]) 
    labels = ms.labels_
    labels_unique = np.unique(labels)
    n_clusters = len(labels_unique)
    print(n_clusters)
    df_clust = pd.DataFrame({"coords":dp_df.iloc[:,3:lt].columns.tolist(), "label": labels})
    print(df_clust['label'])
    clustdir=workdir+'/input/cluster/'
    if not os.path.exists(clustdir):
        os.mkdir(clustdir)
    for cno in range(0,n_clusters):
        clustname=clustdir+'/cluster_'+str(cno)+'_samples.txt'
        df_clust2= df_clust.loc[df_clust.label == cno]
        print(df_clust2)
        df_clust2['coords'].to_csv(clustname, sep='\t',index=False, header=False)
    df_clust.sort_values(by=['label']).to_csv(clustdir+'/Samples_clust.txt', sep='\t',index=False)

def main():
    indx = []
    dfs_cov=[] #coverage dataframe
    dp_df=pd.DataFrame()
    dp_dfn=pd.DataFrame()
    parser = argparse.ArgumentParser(description='CNVexpo ')
    parser.add_argument('--infile', required=True,help='sample list')
    parser.add_argument('--bamfolder', required=True,help='location of the bam files.')
    parser.add_argument('--workdir', required=True,help='location of the work directory.')
    parser.add_argument('--threads', required=True,help='number of cpus.')
    args = parser.parse_args()
    workdir=args.workdir
    fname=workdir+'/input/'+args.infile
    bamfolder=args.bamfolder

    threads= args.threads
    outdir = workdir+'/output/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if os.path.isfile(fname) == False:
        raise Exception("sample list file does not exist,check location")
    lines = [line.strip() for line in open(fname, 'r')]
    prep_fun(lines,bamfolder,workdir,threads )


if __name__ == "__main__":
    main()
