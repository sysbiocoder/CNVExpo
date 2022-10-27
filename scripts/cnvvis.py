# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 11:46:58 2021

@author: sth036
"""
import pandas as pd
import numpy as np
import time
import tkinter
import tkinter as tk
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
import tkinter.filedialog as filedialog
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from pandastable import Table, TableModel
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm 
from matplotlib.patches import Rectangle
import pybedtools
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import socket
import os.path as op
import os
import sys
import cnvplot as cp
import argparse

class cnvexpo(Frame): 
     
    def __init__(self, parent=None):
        Frame.__init__(self, parent)
        global data_table
        self.parent = parent
        #parent.grid_propagate(True)
        self.active = None
        def export(widget):
            global gen
            if len(self.str2.get())==0 :
                messagebox.showinfo("Missing sampleid","Enter sampleid")
                return
            if combo_box.get() == 'select p-value' :
                messagebox.showinfo("Missing p-value","choose p-value")
                return
            if u_df_in > 1:
                sdf=u_df
            else:
                sdf= gene_df.loc[gene_df.zscore <= 0.05]
            #sdf=sdf[sdf.Gene.isin(gen)]

            try:
            # with block automatically closes file
                with filedialog.asksaveasfile(mode='w', defaultextension=".xlsx") as file:
                    sdf.to_excel(file.name)
            except AttributeError:
                # if user cancels save, filedialog returns None rather than a file object, and the 'with' will raise an error
                print("The user cancelled save")
            print(sdf)
            
        def change_df(event):
            global u_df
            #global u_df_in
            #u_df_in = u_df_in+1
            
            zs= combo_box.get()
            #xs= -3.2
            #match zs:
            #    case =0.05:
            #        xs=-3.2
            #     case=0.001:
            #        xs=-5.2
            zpval={
                '0.05':-1.64,
                '0.001':-3.1,
                '0.01':-2.3,
                '0.0001':-3.7
                }
            xs=zpval.get(zs)
            #print(zs)
            tsample=workdir+"/output/report/%s_depth_cal.txt" %self.str2.get()
            u_df=pd.read_csv(tsample,"\t")
            print("udf",u_df)
            u_df.zscore=u_df.zscore.astype(float)
            u_df=u_df.loc[u_df.zscore <= xs]
            
            self.tbl.model.df=u_df
            self.tbl.redraw()
            #print(df['Gene'])
            #sdf=df[df.Gene.isin(gen)]
            
        def selection(txt, widget):
            global count
            
            global df
            if count%2 !=0:
                widget.config(relief="sunken",bg="red")
            if count%2 ==0:
                unselect(txt,widget)
            color=widget.cget("bg").split(",")
            if  'red' in color:
                if txt not in gen:
                    gen.append(txt)
            if 'medium sea green' in color:
                while txt in gen:
                    gen.remove(txt)
            count +=1
            
        def unselect(txt, widget):
            widget.config(relief="groove",bg="medium sea green")
            return

        def get_sample(sampleid, gene_forp):
            global gene_df
            global exon_df
            Gene=gene_forp
            df= pd.read_csv(workdir+"/output/report/{}_depth_cal.txt".format(sampleid),"\t")
            idx=[]
            #df_tmp2=pd.DataFrame(columns=["chr","start","end","exons","strand","Gene"])
            #indfe= pd.read_csv("input/FH-exons-v2.bed","\t",usecols=[0,1,2,3,5,6],names=["chr","start","end","exons","strand","Gene"])
            df_tmp2=pd.DataFrame(columns=["chr","start","end","Gene"])
            indfe= pd.read_csv("data/ccds.gtf","\t",usecols=[0,1,2,3],names=["chr","start","end","Gene"])
            print(Gene)
            #print(indfe)
            dfe=pd.DataFrame(indfe.loc[indfe['Gene'].isin(Gene)])
            dfe=dfe.reset_index()
            dfe = dfe.astype({"chr":"str","start":"int","end":"int","Gene":"str"})
            dfe.drop('index',  axis=1, inplace=True)
            df["Pos"].astype(int)
            #print("dfe",dfe)
            #print("df",df)
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
            c_df.columns =['chr', 'Pos', 'Pos1', 'Gene', 'DP', 'GT', 'MQ', 'DP4', 'start', 'end', 'proportion','Meancoverage', 'std','expecteddepth', 'copynumber', 'zscore']
            c_df.drop('Pos1',  axis=1, inplace=True)
            
            #df_in= df.loc[df['Gene'].isin([Gene])]
            df_in=df
            df_in.drop('Pos1',  axis=1, inplace=True)
            df_ex=c_df
            #print("df_ex",df_ex)
            #print("df_in",df_in)
            gene_df=df
            exon_df=df_ex


        frame1 = Frame(self.parent,background="white smoke")
        frame1.grid(row=0, column=0, sticky=N+S+E+W,rowspan=3)
        screen_width = frame1.winfo_screenwidth() * 0.18
        screen_height = frame1.winfo_screenheight() * 0.5
        tbl_height = frame1.winfo_screenheight() * 0.2
        hm_height = frame1.winfo_screenheight() * 0.006
        L1 = Label(frame1,text="Select genes",width=10, bg='grey23',fg='white')
        L1.grid(row=0, column=1, sticky=N+S+E+W, padx=5, pady=5)
        #frame1.grid_propagate(False)
        global count
        global txt
        global gen
        gen = []
        count=1
        ##For genes
        print(genes)
        for k in range(len(genes)):
            button = Button(frame1,text=genes[k],width=10,bg='medium sea green',relief=GROOVE)
            button.config(height = 2,command=lambda pb=button, x=genes[k]: selection(x,pb))
            button.grid(row=k+1, column=1,sticky=W, padx=5)
            #print(genes[k])
        L2 = Label(frame1,text="Enter sampleid",width=10, bg='grey23',fg='white')
        L2.grid(row=k+2, column=1, sticky=N+S+E+W, padx=5, pady=5)
        T1 = Entry(frame1, width = 10)
        T1.grid(row=k+3, column=1,sticky=W, padx=5)
        
        
        def cnv_plotx(widget):
            global widg
            global sample_id
            global gene_forp
            #widg= widget
            sample_id=T1.get()
            if len(T1.get())==0 :
                messagebox.showinfo("Missing sampleid","Enter sampleid")
            
            if len(gen)==0 :
                messagebox.showinfo("Missing gene","select gene")
            
            gene_forp=gen
            get_sample(sample_id, gene_forp)
            print("gene_forp", gene_forp)
            cp.call_plotx(sample_id,gene_forp,workdir)
        
        def call_igvx(widget):
            import igv
            global igvsample_id
            
            global igvgene_id
            if len(T1.get())==0 :
                messagebox.showinfo("Missing sampleid","Enter sampleid")
            
            if len(gen)==0 :
                messagebox.showinfo("Missing Gene","select gene")
            
            print("IGV called")
            igvsample_id=T1.get()
            igvgene_id=gen[0] #gL.cget("text")
            print("igvgene",gen[0])
            igv.main_fn(igvsample_id,igvgene_id,bamfolder)  
        ## for plot
        button2=Button(frame1,text="Draw plot",width=10,bg='light steel blue',relief=GROOVE)
        button2.config(height = 2,command=lambda plb=button2, : cnv_plotx(plb))
        button2.grid(row=k+4, column=1,sticky=W, padx=5, pady=10)
        
        button3=Button(frame1,text="View BAM",width=10,bg='light steel blue',relief=GROOVE)
        button3.config(height = 2,command=lambda lb=button3, : call_igvx(lb))
        button3.grid(row=k+5, column=1,sticky=W, padx=5, pady=10)

        ##For heatmap
        global frame2
        frame2 = Frame(self.parent,background="white",height=4)  
        frame2.grid(row=0, column=1)
        frame2.grid_propagate(True)
        pbutton=Button(frame2,text="<<",width=10,bg='light steel blue',relief=GROOVE)
        pbutton.config(command=lambda :get_gene("p"))
        pbutton.grid(row=0, column=0, padx=5, pady=2)
        gL = Label(frame2,text=hmgene,width=10, bg='grey23',fg='white')
        gL.grid(row=0, column=1, padx=5, pady=2) 
        nbutton=Button(frame2,text=">>",width=10,bg='light steel blue',relief=GROOVE)
        nbutton.config(command=lambda :get_gene("n"))
        nbutton.grid(row=0, column=2,sticky=W, padx=5, pady=2)
        
        ##For heatmap
        global frame3
        frame3 = Frame(self.parent,background="white smoke", width=10, height=4)  
        frame3.grid(row=2, column=1, sticky=W)
        frame1.grid_propagate(True)
        scroll_y = tk.Scrollbar(frame3, orient="vertical")
        scroll_y.pack(side="right", expand=True, fill="y")
        
        ##For table
        frame4 = Frame(self.parent,background="white")  
        frame4.grid(row=3 ,column=0, sticky=N+S+E+W, columnspan=6)
        firstsample=workdir+"/output/report/%s_depth_cal.txt" %samples[0]
        gene_df=pd.read_csv(firstsample,"\t")
        #gene_df.columns=["chr","Pos","DP","GT","MQ","Allele fraction","start","end","Gene","expecteddepth","copynumber","Zscore"]
        self.tbl=gtb = Table(frame4, dataframe=gene_df,height=tbl_height, width = screen_width)
        self.tbl.grid(row=1, column=2,rowspan=8)
        self.tbl.grid_propagate(True)
        gtb.show()
        
        ## For download
        def refresh_df(event):
            combo_box.current(0)
            #print("djdk")
        frame5 = Frame(self.parent,background="white")
        frame5.grid(row=4, column=1, columnspan=2)
        frame5.grid_rowconfigure(0, weight=1)
        
        tL2 = Label(frame5,text="Enter sampleid",width=15, bg='grey23',fg='white')
        tL2.grid(row=0, column=1, sticky=W, padx=5, pady=5)       
        self.str2 =StringVar()
        
        self.e2 = tk.Entry(frame5, textvariable=self.str2)
        #self.e2.bind('<Tab>', (lambda event: on_changed))
        self.e2.grid(row=0, column=2, sticky=W, padx=5, pady=5)
        
        self.str2.trace("w", lambda name, index,mode, var=self.str2: refresh_df(self.str2))
        #self.e2.bind_all('<<Modified>>', lambda event: on_changed)

        ##combobox to filter
        L3 = Label(frame5,text="Filter by p-value",width=15, bg='grey26',fg='white')
        L3.grid(row=0, column=3, sticky=W, padx=5, pady=5)
        combo_choices = ['select p-value','0.05','0.01', '0.001','0.0001']
        choice = StringVar()
        combo_box = ttk.Combobox(frame5, textvariable="select p-value")
        combo_box['values'] = combo_choices
        combo_box.current(0)
        combo_box.grid(column=4, row=0, sticky=W,padx=10)
        print(gen)
        combo_box.bind('<<ComboboxSelected>>',change_df) 

        button4=Button(frame5,text="Download table",width=15,bg='light steel blue',relief=GROOVE)
        button4.config(height = 2,command=lambda mb=button4, : export(mb))
        button4.grid(row=0, column=5)


        def checkcanvas(pcounter):
            try:
                hcanvas.get_tk_widget().pack_forget()
                pcounter=pcounter +1
                print("canvas  exist")
            except:
                #get_gene('b')
                pass
        def draw_hmap(pcounter):
            global figh
            global ax_h
            global hcanvas, hsamples
            checkcanvas(pcounter)
            figh,ax_h = plt.subplots(2,figsize=(8,hm_height))
            hcanvas = FigureCanvasTkAgg(figh, master= frame3)
            hcanvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.X, expand=0)
            hsample=samples[0]
            hsamples=samples
            sam_df=pd.read_csv(workdir+"/output/report/{}_depth_cal.txt".format(hsample),"\t")
            sam_df.columns=["chr","Pos","Gene","DP","GT","MQ","Allele fraction","start","end","proportion","Meancoverage","std","expecteddepth","copynumber","zscore"]
            sam_df['zscore']=-sam_df['zscore'].abs()
            print("sam_df",sam_df)
            gen_df=sam_df.loc[sam_df['Gene'].isin([hmgene])]
            Position=gen_df['chr'].values[0]
            cnv_df=gen_df[['Pos','copynumber']]
            cnv_df.columns=[Position,hsample]
            plot_df=gen_df[['Pos','zscore']]
            plot_df.columns=[Position,hsample]
            hm=2
            for hi in range(1,len(hsamples)):
                get_df= pd.read_csv(workdir+"/output/report/{}_depth_cal.txt".format(hsamples[hi]),"\t")
                cnv_df.insert(hm,hsamples[hi],get_df[["copynumber"]]) #copynumer
                plot_df.insert(hm,hsamples[hi],get_df[["zscore"]])  #zscore
                hm=hm+1
            plot_df = plot_df.dropna(axis=0)
            cnv_df = cnv_df.dropna(axis=0)
            plot_df.set_index(Position, inplace=True, drop=True)
            cnv_df.set_index(Position, inplace=True, drop=True)
            plot_df=plot_df.astype(float)
            cnv_df=cnv_df.astype(float)
    
            # Creating plot
            hn=len(hsamples)
            colormap=sns.color_palette("mako", 2,as_cmap=True)
            sns.set(font_scale=0.6)
            bounds = [-5, -3.1,-2.3,-1.64, 0]
            norm=BoundaryNorm(bounds, colormap.N) 
            res = sns.heatmap(plot_df.T,cmap="cividis",ax=ax_h[0], xticklabels= False,center=-3.1, vmin=-4, vmax=4,norm=norm)
            res.tick_params( bottom=False)
            res.set_xlabel('')
            res.set_title('Z-score')
            res.set_yticklabels(res.get_yticklabels(),rotation=0)
            colormap1=sns.color_palette("bwr_r",3,as_cmap=True)
            bounds2 = [2.4, 2, 1.6]
            norm2=BoundaryNorm(bounds2, colormap.N) 
            res1 = sns.heatmap(cnv_df.T,cmap=colormap1, ax=ax_h[1],center=2, vmin=0, vmax=4, robust=True)

            res1.set_yticklabels(res1.get_yticklabels(),rotation=0)
            res1.set_title('CNV')
                # show plot
            figh.tight_layout()
            #get_gene('b')
        

        def get_gene(prev):
            global hmgene
            global pcounter
            pcounter=0
            
            for k in range(len(genes)):
                if(genes[k]) == gL.cget("text") :
                    gcounter= k
                    
            if(prev=="b"):
                hmgene=genes[0]   
                print(hmgene)
                gL["text"]=genes[0]
                pcounter=pcounter+1
                #draw_hmap(pcounter)
            if((prev=="p") & (gcounter >= 1)):
                gcounter=gcounter-1
                hmgene=genes[gcounter]
                gL["text"]=genes[gcounter]
                print(hmgene)
                pcounter=pcounter+1
                #draw_hmap(pcounter)

                
            if((prev=="n") & (gcounter < len(genes)-1)) :
                #print(pcounter)
                #print(gcounter)
                #print(len(genes))
                gcounter=gcounter+1
                hmgene=genes[gcounter]
                gL["text"]=genes[gcounter]
                print(hmgene)
                pcounter=pcounter+1
                
            if((prev=="n") & (gcounter >= len(genes))) :
                gcounter=0
                hmgene=genes[gcounter]
                gL["text"]=genes[gcounter]
                print(hmgene)
                pcounter=0
            draw_hmap(pcounter)
        draw_hmap(0)
              

#Create dataframe
root = tkinter.Tk()
root.title("CNV Expo")
parser = argparse.ArgumentParser(description='CNVexpo ')
parser.add_argument('--genelist', required=True,help='gene list')
parser.add_argument('--samplelist', required=True,help='sample list (maximum 5 samples allowed')
parser.add_argument('--workdir', required=True,help='location of the work directory.')
parser.add_argument('--bamfolder', required=True,help='location of the bam directory.')
args = parser.parse_args()
workdir=args.workdir
gname=workdir+'/input/'+args.genelist
sname=workdir+'/input/'+args.samplelist
global bamfolder
bamfolder=args.bamfolder

with open(gname) as ginp:
        genes = []
        for gline in ginp:
            genes.append(gline.rstrip('\n'))
hmgene=genes[0]

with open(sname) as sinp:
        samples = []
        for sline in sinp:
            samples.append(sline.rstrip('\n'))
count =1
u_df_in=1


root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
root.maxsize(1500,10000)
root.resizable(True,True)
root.configure(background="white")
display = cnvexpo(root)
root.mainloop()

