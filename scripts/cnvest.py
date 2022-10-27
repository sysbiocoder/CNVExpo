import os
import pysam
import pysam.bcftools as bcftools
import numpy as np
import pandas as pd
from cyvcf2 import VCF
import warnings
warnings.filterwarnings('ignore')
def cnv_cal(sam_lst,dfs_cov,indx):
    #print("dfs_cov",dfs_cov)
    df_sam = pd.DataFrame({'idx':indx, 'dfs':dfs_cov}) 
    #print("df_sam",df_sam['dfs'].iloc[0][["chr"]])
    dp_df = pd.DataFrame()  #total sample dataframe
    dp_df["chr"]=df_sam['dfs'].iloc[0]["chr"]
    dp_df["start"]=df_sam['dfs'].iloc[0][["start"]]
    dp_df["end"]=df_sam['dfs'].iloc[0][["end"]]
    dp_df["Pos"]=df_sam['dfs'].iloc[0][["Pos"]]
    dp_df["Gene"]=df_sam['dfs'].iloc[0][["Gene"]]
    #Proportion table
    dp_dfn=pd.DataFrame()
    for dpidx in range(0,len(indx)):
        sam_name=sam_lst[dpidx]
        dp_dfn=pd.DataFrame({'chr':df_sam['dfs'].iloc[dpidx]["chr"],'Pos':df_sam['dfs'].iloc[dpidx]["Pos"],'DP':df_sam['dfs'].iloc[dpidx]["DP"]})
        dp_df=pd.merge(dp_df,dp_dfn,on=['chr','Pos'],how='outer')
        dp_df.rename(columns={ dp_df.columns[dpidx+5]: sam_lst[dpidx] }, inplace = True)
        dp_df.drop_duplicates(subset=['chr', 'Pos'],keep='first',inplace=True)
        dp_df.replace([np.inf, -np.inf], 0, inplace=True)
        dp_df.replace([np.nan, -np.nan], 0, inplace=True)
        dp_df.fillna(0)
        dp_df=dp_df.reset_index(drop=True)
    locn= 5 + len(indx)  #chr,start,end,pos,Gene+no of samples
    locn2= locn + len(indx)
    posn= 5+dpidx #chr,start,end,pos,Gene+sampleid  
    prop_df=pd.DataFrame(dp_df)
    headers=pd.DataFrame()
    #headers =  prop_df.iloc[:,5:locn].groupby(prop_df.iloc[:, 4]).sum()
    headers1 =  pd.DataFrame(prop_df.iloc[:,5:locn].sum())
    #headers = headers.rename_axis('Gene').reset_index()
    #prop_df2 = prop_df.merge(headers,on=['Gene'])
    #print("prop_df2",prop_df2) #595781
    #prop_df_div = prop_df2.iloc[:,np.r_[0:5,locn:locn2]]
    prop_df_div = headers1

    #print("prop_df",prop_df.iloc[:,5:locn]) #595781
    #print("prop_df_div",prop_df_div) #595781

    #prop_df_div.columns=prop_df.columns
    #prop_df.iloc[:,5:locn]=prop_df.iloc[:,5:locn].astype(int).divide(prop_df_div.iloc[:,5:locn].astype(int),axis=1)
    prop_df.iloc[:,5:locn]=prop_df.iloc[:,5:locn].astype(int).div(prop_df_div.iloc[:,0],axis=1)
    prop_df['Meancoverage']= prop_df.iloc[:,5:locn].mean(numeric_only=True,axis=1)#/float(len(indx)) #sum of each position
    
    prop_df['std']=prop_df.iloc[:,5:locn].std(numeric_only=True, axis=1)
    dp_df=pd.merge(dp_df,prop_df,on=['chr','Pos','start','end','Gene'])
    #print("prop_df",prop_df)

    #CNV calculation
    for cpidx in range(0,len(indx)):
        sam_name=sam_lst[cpidx]
        df_sam['dfs'].iloc[cpidx].drop_duplicates(subset=['chr', 'Pos'],keep='first',inplace=True)
        dindx=5+cpidx
        lindx=5+len(indx)
        dp_df=prop_df.iloc[:,[0,1,2,3,4,dindx,lindx,lindx+1]]
        sam_df = pd.DataFrame(df_sam['dfs'].iloc[cpidx]) 
        df_sam['dfs'].iloc[cpidx]=pd.merge(sam_df,dp_df,on=['chr','Pos','start','end','Gene'])
        df_sam['dfs'].iloc[cpidx]['expecteddepth']= df_sam['dfs'].iloc[cpidx]["DP"].div(df_sam['dfs'].iloc[cpidx][sam_name])*(df_sam['dfs'].iloc[cpidx]["Meancoverage"])
        df_sam['dfs'].iloc[cpidx]["copynumber"]=2*(df_sam['dfs'].iloc[cpidx][sam_name].astype(float)).div(df_sam['dfs'].iloc[cpidx]["Meancoverage"].astype(float))
        df_sam['dfs'].iloc[cpidx]["zscore"]= ((df_sam['dfs'].iloc[cpidx][sam_name]).sub(df_sam['dfs'].iloc[cpidx]["Meancoverage"])).div(df_sam['dfs'].iloc[cpidx]['std'])
        
        df_sam['dfs'].iloc[cpidx].drop(df_sam['dfs'].iloc[cpidx][df_sam['dfs'].iloc[cpidx]['expecteddepth'] < 20].index,inplace=True)
        #print(sam_name,df_sam['dfs'].iloc[cpidx])

    return(df_sam['dfs'],dp_df)