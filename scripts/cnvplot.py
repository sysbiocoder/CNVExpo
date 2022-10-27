# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 11:40:02 2021

@author: sth036
"""
import pandas as pd
import numpy as np
import math
import tkinter
import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.patches import Rectangle
import numpy.ma as ma
import pybedtools


def get_inputsamplex(sampleid,workdir):
    global gene_df
    #global exon_df
    global df_ex
    global df_in
    
    Gene=Genex #"LDLR"
    #sampleid="22143"
    df= pd.read_csv(workdir+"/output/report/{}_depth_cal.txt".format(sampleid),"\t")
    idx=[]
    df_tmp2=pd.DataFrame(columns=["chr","start","end","Gene"])
    indfe= pd.read_csv("data/ccds.gtf","\t",usecols=[0,1,2,3],names=["chr","start","end","Gene"])
    dfe=pd.DataFrame(indfe.loc[indfe['Gene'].isin([Gene])])
    dfe=dfe.reset_index()
    dfe = dfe.astype({"chr":"str","start":"int","end":"int","Gene":"str"})
    dfe.drop('index',  axis=1, inplace=True)
    df["Pos"].astype(int)
    df=pd.DataFrame(df.loc[df['Gene'].isin([Gene])])
    #print("df",df)
    #print("dfe",dfe)
    indx=1
    row=1
    #print("loc")
    #print(dfe['start'].iloc[indx],dfe['end'].iloc[indx])
    #print(df['Pos'].iloc[row])
    df_ex = pd.DataFrame()
    #print(len(dfe))
    #print(len(df))
    df.insert(2,'Pos1',df['Pos'].astype(int) + 1)
    #print(df.head)
    a=pybedtools.BedTool.from_dataframe(dfe)
    b=pybedtools.BedTool.from_dataframe(df)
    #print("ab",a,b)
    c= a.intersect(b,wb=True)
    
    c_df=pybedtools.BedTool.to_dataframe(c,disable_auto_names=True, header=None)
    c_df = c_df.iloc[:,4:]
    #print("c_df",c_df)
    c_df.columns =["chr","Pos","Pos1","Gene","DP","GT","MQ","DP4","start","end","proportion","Meancoverage","std","expecteddepth","copynumber","zscore"]
    c_df.drop('Pos1',  axis=1, inplace=True)
    #print("c_df",c_df)
    df_in= df.loc[df['Gene'].isin([Gene])]
    df_in.drop('Pos1',  axis=1, inplace=True)
    df_ex=c_df
    df_ex.dropna(inplace=True)
    df_in.dropna(inplace=True)
    #print("exons",df_ex)
    #print("genes",df_in)

    #print(list(df_ex.columns))
    return(df_ex,df_in)
def get_exonlocx(Genex):
    global df2
    df3= pd.read_csv("data/ccds.gtf","\t",header=None)
    df3.columns=["chr","start","end","Gene"]
    df2=df3.loc[df3['Gene'].isin([Genex])]
    #print(df2)
    return(df2)


def get_cnvlocx(Genex):
    #print("cnvlic", Genex)
    global df1,dfi1
    df1=df_ex#.loc[df_ex['Gene'].isin([Genex])]
    #print(df1['Pos'])
    dfi1=df_in#.loc[df_in['Gene'].isin([Genex])]
     #print(dfi1['pos'])
    return(df1,dfi1)

def handler1(Genex):
    print("Gene.",Genex)
    df2=get_exonlocx(Genex)
    return(df2)
           
def handler2(Genex,Samplex,workdir):
    #Samplex="22143"
    #print("Sample.",Samplex)
    (dfs_ex,dfs_in)=get_inputsamplex(Samplex,workdir)
    #print("callcnvlocx")
    (dfc_ex,dfc_in)=get_cnvlocx(Genex)
    return(dfc_ex,dfc_in)
    
def draw_plot2(ndf,getex,Samplex):
    root2 = tk.Toplevel()
    fig,ax = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [4,1]})
    pcanvas = FigureCanvasTkAgg(fig, root2)
    #canvas.show()
    pcanvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    toolbar = NavigationToolbar2Tk(pcanvas, root2)
    toolbar.update()
    pcanvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    scrollbar = tkinter.Scrollbar(root2, orient=tk.HORIZONTAL)
    scrollbar.pack(side=tk.BOTTOM, fill=tk.X)
    scrollbar["command"] = pcanvas.get_tk_widget().xview
    pcanvas.get_tk_widget()["xscrollcommand"] = scrollbar.set
            
    #canvas._tkcanvas.pack(anchor=tk.W, side=tk.TOP, fill=tk.BOTH, expand=1)
    #Samplex="22143"
    #plt.tight_layout()
    
    ndf.rename( columns={'zscore' :'Zscore'}, inplace=True )
    print(ndf)
    ndf['copynumber']=pd.to_numeric(ndf['copynumber'], errors='coerce') 
    ndf['Zscore']=pd.to_numeric(ndf['Zscore'], errors='coerce') 
    px=ndf.iloc[:,1]  #posn
    y1=ndf.iloc[:,13]  #cN
    y2=ndf.iloc[:,14]     #zscore
    #print("CN",px)
    print(len(y1))
    y1a = ma.array(y1)
    y2a = y1.copy()
    za=ma.array(y2)
    #print("y1a",y1a)
    y3a = y1.copy()
    ax[0].plot(px, y1a, color="grey", marker="o")
    deln_threshold = 1.6
    zs_threshold = -1.64
    deln = np.ma.masked_where(y1a > deln_threshold, y1a)
    ax[0].plot(px, deln, color="red", marker="o")
    
    dup_threshold = 2.4
    dup = np.ma.masked_where(((y2a < deln_threshold) | (y2a < dup_threshold) ) , y2a)
    ax[0].plot(px, dup, color="blue", marker="o")
    ax[0].axhline(y=1.6, color='r', linestyle='dotted')
    ax[0].axhline(y=2.0, color='g', linestyle='--')
    ax[0].axhline(y=2.4, color='b', linestyle='dotted')
    # set x-axis label
    s= ''.join(list(set("hg38 "+df1["chr"])))
    ax[0].set_ylim([0,11])
    # set y-axis label
    ax[0].set_title('Sample= '+str(Samplex)+', Gene = '+str(Genex)+', Region = '+str(getex))
    ax[0].set_ylabel("Copynumber",color="red",fontsize=14)
    
    #ey=ndf.iloc[:,7]#expected depth
    #ley=ey.apply(lambda x: np.log2(x))#np.log2(ey)
    #ax[0].plot(px,ley , color="grey", linestyle="--")
    
    ax2=ax[0].twinx()
    y2b=ma.array(y2)
    pvalm=np.ma.masked_where(-abs(y2b) <= -2.3 , -abs(y2b))
    ax2.plot(px, -abs(y2b),color="purple",marker="+")
    ax2.plot(px, -abs(pvalm),color="black",marker="+")
    ax2.axhline(y=-2.3, color='grey', linestyle='dotted')
    ax2.axhline(y=-1.64, color='grey', linestyle='--')

    ax2.set_ylabel("-abs(zscore)",color="blue",fontsize=14)
    ax2.set_ylim([-10,1])
    

    
    py=y1
    line1=ax[1].plot(px,py,color="white")
    ln = line1.pop(0)
    ln.remove()
    ax[1].set_yticks((0,1))
    ax[1].set_yticks([])
    ax[1].set_xlabel(s,fontsize=8)
    xmin=min(px)
    xmax=max(px)
    count_row = df2.shape[0]
    #print(count_row)
    for pk in range(0,count_row):
        m=df2.iloc[pk,1]
        n=df2.iloc[pk,2]
        exon=n-m+1
        p=df2.iloc[pk,2]
        j=pk+1
        if(j < count_row):
            q=df2.iloc[j,1]
            intron= q - p
        else:
            xmid= p + intron/2
            q=xmid
            
        xmid= p + intron/2
        ymid=2
        xval=[p,xmid,q]
        yval=[1,ymid,1]
        ax[1].plot(xval,yval,color="grey")
        ax[1].add_patch(Rectangle((m,0),exon,3,color="grey"))
                
                #ax.plot([1,5,2],[2,3,4],color="cyan")
                #ax.add_patch(Rectangle((2, 2), 1, 3,color="yellow"))
        
    pcanvas.draw()

def call_plotx(sample_id,gene_listx,workdir):
    global Genex
    #gene_listx=["LDLR"]
    for pi in range(len(gene_listx)):
        #cnvtab4.widg.config(relief="groove",bg="grey23")
        Genex=gene_listx[pi]
        df1=handler1(Genex)
        
        (df1,df1i)=handler2(Genex,sample_id,workdir)
        
        draw_plot2(df1,"exon",sample_id) 
        draw_plot2(df1i,"gene",sample_id)

            


