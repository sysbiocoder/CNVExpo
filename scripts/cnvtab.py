import pandas as pd
import numpy as np
import scipy.stats as sp
import math 

def tab_form(df_tab,sample,workdir):
    df_tab['MQ'] = pd.to_numeric(df_tab['MQ'])
    df_tabx = pd.DataFrame()
    df_tabh=pd.DataFrame(columns=["chrno","start","end","Gene","copynumber","Quality","length","genotype","combined pvalue","decipher"])
    
    genes=df_tab['Gene'].unique()
    for gns in range(0,len(genes)):
        df_tabx=df_tab[df_tab['Gene'] ==genes[gns]]
        df_tabx=df_tabx.reset_index()
        
        df_tabx['diff_pos'] = df_tabx.Pos - df_tabx.Pos.shift()
        #print("pos",df_tabx.Pos)
        #print("shift",df_tabx.Pos.shift().dropna().astype(int))
        diff= df_tabx['diff_pos'].mean()#quantile(0.2)
        #print(diff)
        df_tabx['diff2']=df_tabx['copynumber'].diff() < 0
        diff2=df_tabx[df_tabx.diff2=='False']
        #print(diff2.index)
        #df_tabx['bool_eq']=df_tabx.diff_pos.lt(int(diff))
        #print(df_tabx['bool_eq'])
        #bins1 = df_tabx[df_tabx.bool_eq ==0]#.tolist()
        #index=bins1.index
        df_tabx['bool_eq']=df_tabx.diff_pos.lt(diff)
        bins1 = df_tabx[df_tabx.bool_eq ==0]
        index=bins1.index
        #print(index)
        bin2=[]
        k=1
        j=1
        for i in range(0,len(index)-1):
            j=index[i]
            k=index[i+1]-1
            
            if df_tabx.loc[j]['Pos'] != df_tabx.loc[k]['Pos']:
                bin2.append(df_tabx.loc[j]['Pos'])
                bin2.append(df_tabx.loc[k]['Pos'])
    

        df_tab2 = pd.DataFrame()
        df_tab2['chrno']=df_tabx.groupby(pd.cut(df_tabx['Pos'],bin2))['chr'].first()
        df_tab2['start']=df_tabx.groupby(pd.cut(df_tabx['Pos'],bin2))['Pos'].min()
        df_tab2['end']=df_tabx.groupby(pd.cut(df_tabx['Pos'],bin2))['Pos'].max()
        gname=df_tabx['Gene'].unique()[0]
        df_tab2['Gene']=str(gname)
        cn=df_tabx.groupby(pd.cut(df_tabx['Pos'],bin2))['copynumber'].mean()
        df_tab2['copynumber']=cn.to_frame().values
        df_tab2['Quality']=df_tabx.groupby(pd.cut(df_tabx['Pos'],bin2))['MQ'].mean().to_frame()
        pos1=df_tabx.groupby(pd.cut(df_tabx['Pos'],bin2))['Pos'].min().to_frame()
        pos2=df_tabx.groupby(pd.cut(df_tabx['Pos'],bin2))['Pos'].max().to_frame()
        df_tab2['length']=df_tab2['end']-df_tab2['start']
        lent=k-j
        
        df_tab2['genotype']=df_tabx.groupby(pd.cut(df_tabx['Pos'],bin2))['GT'].unique()
        df_tab2['combined pvalue']=df_tabx.groupby(pd.cut(df_tabx['Pos'],bin2))['zscore'].apply(lambda c: c.abs().sum())/math.sqrt(abs(lent))
        df_tab2=df_tab2.iloc[::2]
        df_tab2.reset_index(drop=True, inplace=True)
        df_tab2['copynumber'] = pd.DataFrame(df_tab2['copynumber'])
        region=[]
        af=[]
        afpos=[]
        df_tab2=df_tab2.replace('chr', '',regex=True)
        pd.set_option("display.precision", 2)
        dlink='https://www.deciphergenomics.org/browser#q/'
        for index in df_tab2.index:
            if(df_tab2['copynumber'][index] <= 1.6):
                region.append("<a href="+'\''+dlink +str(df_tab2['chrno'][index])+':'+str(df_tab2['start'][index])+'-'+str(df_tab2['end'][index])+'/location/'+str(df_tab2['chrno'][index])+':'+str(df_tab2['start'][index])+'-'+str(df_tab2['end'][index])+'\''+'>'+ 'g'+'.'+str(df_tab2['start'][index])+'_'+str(df_tab2['end'][index])+'del'+'</a>' )
            if(df_tab2['copynumber'][index] >= 2.4):
                region.append("<a href="+'\''+dlink +str(df_tab2['chrno'][index])+':'+str(df_tab2['start'][index])+'-'+str(df_tab2['end'][index])+'/location/'+str(df_tab2['chrno'][index])+':'+str(df_tab2['start'][index])+'-'+str(df_tab2['end'][index])+'\''+'>'+ 'g'+'.'+str(df_tab2['start'][index])+'_'+str(df_tab2['end'][index])+'dup'+'</a>' )
            df_tab2['combined pvalue'][index]=pd.Series(sp.norm.sf(df_tab2['combined pvalue'][index]/np.sqrt(df_tab2['length'][index])))
        #print(sample)
        #print(df_tab2)
        df_tab2.reset_index(drop=True, inplace=True)
        df_tab2['decipher']=pd.Series(region)
        df_tab2.drop(df_tab2[df_tab2['length'] < 50].index, inplace = True)
        df_tabh=pd.concat([df_tabh,df_tab2])
        df_tabh.reset_index(drop=True, inplace=True)
        
        
    return(df_tabh)


def make_tab(df_tab,sample,workdir):
    
    df_tab.dropna( inplace = True)
    #print("df_tab",df_tab)
    df_tab.drop(df_tab[df_tab['zscore'].abs().astype(float) < 1.65].index, inplace = True)

    df_tab.drop(df_tab[((df_tab['copynumber'] < 2.4) & (df_tab['copynumber'] > 1.6))].index, inplace = True)
    #df_tab.drop(df_tab[df_tab['MQ'].astype(int) < 60].index, inplace = True)
    df_table1=df_tab[df_tab['copynumber'] > 2.4]
    df_table2=df_tab[df_tab['copynumber'] < 1.6]
    
    df_tabh1=tab_form(df_table1,sample,workdir)
    df_tabh2=tab_form(df_table2,sample,workdir)
    df_tabh=pd.concat([df_tabh1,df_tabh2],axis=0)
    print(sample,df_tabh)
    htmlfile=workdir+'/output/report/'+str(sample)+'output.depth.html'
    df_tabh.to_html(htmlfile, index=False,render_links=True,escape=False)
