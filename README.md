# CNVExpo <br>
``CNVExpo: 
Python package for detection of copy number variation analysis in next generation sequencing data 

CNVExpo uses mapped bam file as input and determine read-depth in each position of the panel or exome bed file provided by user to compute 
CNV. Bcftools is used to determine the read depth per position. Clustering is done based on a control gene, where CNV is not present. Here is the structure of files in CNVExpo.


├── cnvcal.py <br>
├── data <br>
│   ├── ccds.gtf <br>
│   ├── control.bed <br>
│   ├── GRCh38_full_analysis_set_plus_decoy_hla.fa  (GRCh38_full_analysis_set_plus_decoy_hla.fa (ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa) <br>
│   └── GRCh38_full_analysis_set_plus_decoy_hla.fa.fai <br>
└── scripts <br>
    ├── cnvest.py <br>
    ├── cnvplot.py <br>
    ├── cnvtab.py <br>
    ├── cnvvis.py <br>
    ├── findclust.py <br>
    ├── igv.py <br>
    ├── workdir <br>
        ├── input <br>
            ├── bed file <br>
            ├── gene list <br>
            ├── sample list <br>
        ├── output <br>

## Installation <br>
CNVExpo requires the following libraries to be installed and is checked in Linux environment with the version provided below<br>
Requirements: <br>
Python==3.7.13<br>
pysam==0.19.0<br>
pandas==1.3.5<br>
bcftools==1.9 <br>
cyvcf2==0.30.15 <br>
click==7.1.2 <br>
numpy==1.21.5 <br>
scipy==1.7.3 <br>
sklearn==0.0 <br>
vcf-parser==1.6 <br>
pybedtools==0.8.1<br>

The following files are required to be kept in folder data GRCh38_full_analysis_set_plus_decoy_hla.fa, GRCh38_full_analysis_set_plus_decoy_hla.fa.fai 

To install use command <br>

git clone https://github.com/sysbiocoder/CNVExpo.git <br>

Run from folder CNVExpo  <br>
 

##  Step 1: Input requirements<br>
Create a project_folder, which would be the working directory<br>
Inside project folder make a directory input and store bedfiles, sample list and gene list.<br>
CNVExpo also requires location of the bam file.<br>
•	Generating genes/control bedfile:<br>
The bed file is supposed to have 5 fields for each gene or regions delimited by <TAB> space.<br>
The fields required in the bedfile are given below:<br>
Chr start end INFO/Genename strand biotype <br>
<br>
Example bed file (control.bed): <br>
chr1	11107485	11107500	MTOR	-	protein_coding<br>
chr1	11108180	11108286	MTOR	-	protein_coding<br>
chr1	11109289	11109370	MTOR	-	protein_coding<br>
chr1	11109648	11109729	MTOR	-	protein_coding<br>
<br>
•	Generating sample list <br>
Make a sample list text file with each sample per line<br>
<br>
Example sample list (sample.txt)<br>
22101<br>
22102<br>
22105<br>
<br>
•	Generating gene list<br>
Make a gene list text file with each gene per line<br>
<br>
Example gene list (gene.txt)<br>
      PCSK9<br>
      LDLR<br>
      APOB<br>
<br>      
## Step 2: Clustering analysis<br>
To determine the background dataset, findclust.py script could be utilized.<br>
It requires control bed file, bamfolder- the location of bam files, sample list, working directory.<br>
Run the python script as below<br><br>
*python scripts/findclust.py --infile samplet.txt --workdir project_folder  --bedfile test.bed --bamfolder bamfiles_locn –-threads number*<br><br>
<br>
It generates cluster folder inside input directory and different cluster lists with samples per each line in the cluster, which could be used as input for the cnv estimation step.<br>
<br>
<br>
##  Step 3: CNV Estimation<br>
To estimate CNV, cnvcal.py script is utilized.<br>
It requires bed file of the genes/region, bamfolder- the location of bam files, sample list (or the cluster file from previous step copied to the input directory), working directory.<br>
Run the python script as below<br><br>
*python cnvcal.py --infile clusterfile --bedfile target.bed --bamfolder bamfiles_locn --workdir project_folder –-threads number*  <br><br>
It generates vcf files, html report  and depth_cal.txt for each samples.<br>
<br>
<br>
##  Step 4: CNV Visualization<br>
To visualize and explore CNV, cnvvis.py can be used<br>
It requires gene lists, sample lists (Maximum 10 samples) bamfolder- the location of bam files, and the working directory.<br>
Run the python script as below<br><br>
*python cnvvis.py --genelist genes_list.txt --samplelist cluster_1_samples --bamfolder bamfiles_locn --workdir project_folder*<br><br>
It generates visualization tool for the gene list and sample list provided. <br>

