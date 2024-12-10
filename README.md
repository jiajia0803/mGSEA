# mGSEA: metagenomic single sample Gene Set Enrichment Analysis  

mGSEA is a novel tool developed based on the principle of ssGSEA.
mGSEA is a tool that can be applied to the analysis of metagenomic or metatranscriptomic data.
mGSEA can calculate the enrichment scores of BGCs in individual metagenomic samples.
For mGSEA analysis, two key files are required: the abundance data of BGCs and the classification information of BGCs.

# Install deepBGC using conda.

```bash
#Install the deepBGC tool (for predicting BGCs in the gut microbiome. If there is a better tool available, replace this one).
#Add the bioconda channel and conda-forge.
conda config --add channels bioconda
conda config --add channels conda-forge
#Create a virtual environment for deepBGC and install dependent tools.
conda create -n deepbgc python=3.7 hmmer prodigal
#Activate the deepBGC virtual environment.
conda activate deepbgc
#Install deepBGC.
conda install deepbgc
#Download the pre-trained models and Pfam database.
deepbgc download
#Check the downloaded dependencies and models.
deepbgc info
```
# Predict BGCs with the deepBGC tool and filter the prediction results

```bash
#The deepBGC tool is used to predict BGCs. "file_scaffolds.fasta" is the contig file after assembly. (Adding "--prodigal-meta-mode" can detect genes in short contigs).
deepbgc pipeline --prodigal-meta-mode  file_scaffolds.fasta
For each sample, a result file with the suffix xx.bgc.tsv will be generated after deepBGC prediction. To avoid confusion in batch processing, add the sample name to the first column of each xx.bgc.tsv file. "^" indicates the first column, and "sed -i" means to modify the file directly without printing it out.
for file in xx.bgc.tsv; do sed -i "s/^/$file/" "$file"; done
Merge all deepBGC prediction result files with the suffix .bgc.tsv for batch processing, and generate a summary file named all.bgc.tsv.
find ./ -name "*.bgc.tsv"| xargs cat > all.bgc.tsv
#Delete rows with empty values in the 18th column (product_class column) of the all.bgc.tsv file, and output the new results to the all.bgc_1.tsv file.
awk -F, '$11 != ""' all.bgc.tsv > all.bgc_1.tsv
#Filter the results in the all.bgc_1.tsv file where the 12th column (deepbgc_score column) has a value greater than 0.5, and output to the all.bgc_2.tsv file. The all.bgc_2.tsv file contains the BGC prediction results from deepBGC, including the names of the BGC gene sets and possible predicted classifications.
awk -F, '$3 > 0.5' all.bgc_1.tsv > all.bgc_2.tsv
```
# Construction of deepBGC Gene Sets

```bash
#Combine all the genbank files from the deepBGC prediction results into a new folder named "all_deepBGC_gbk" for batch processing.
find -name '*.gbk' -exec cp {} /mnt/data1/ZhouJiaJia/HMP_IBD/HMP_1338/MG_1338/all_deepBGC_gbk \;
#Check which genbank files are in the all_deepBGC_gbk folder.
cd /mnt/data1/ZhouJiaJia/HMP_IBD/HMP_1338/MG_1338/all_deepBGC_gbk && ls
# Use the pick_fa.py script to extract nucleotide sequences from the genbank files.
python pick_fa.py
# The following are some examples.
```
#python
#Write a Python script named pick_fa.py to extract nucleic acid sequences from genbank files such as "SRR5935758_all.1.bgc.gbk", "SRR5935758_all.2.bgc.gbk", "SRR5946523_all.1.bgc.gbk"... genbank files.
from Bio import SeqIO
gbk_filename = "SRR5935758_all.1.bgc.gbk"
faa_filename = "SRR5935758_all.1.bgc.fna"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    output_handle.write(">%s %s\n%s\n" % (
           seq_record.id,
           seq_record.description,
           seq_record.seq))

output_handle.close()
input_handle.close()
gbk_filename = "SRR5935758_all.2.bgc.gbk"
faa_filename = "SRR5935758_all.2.bgc.fna"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    output_handle.write(">%s %s\n%s\n" % (
           seq_record.id,
           seq_record.description,
           seq_record.seq))

output_handle.close()
input_handle.close()
gbk_filename = "SRR5946523_all.1.bgc.gbk"
faa_filename = "SRR5946523_all.1.bgc.fna"
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    output_handle.write(">%s %s\n%s\n" % (
           seq_record.id,
           seq_record.description,
           seq_record.seq))

output_handle.close()
input_handle.close()
```bash
#Combine all the nucleotide sequence files extracted from the genbank files into a new folder named "all_deepBGC_fna" for batch processing.
find -name '*.bgc.fna' -exec cp {} /mnt/data1/ZhouJiaJia/HMP_IBD/HMP_1338/MG_1338/all_deepBGC_fna \;

#Check which nucleotide sequence files are in the all_deepBGC_fna folder.
cd /mnt/data1/ZhouJiaJia/HMP_IBD/HMP_1338/MG_1338/all_deepBGC_fna && ls

#Add the filename to the beginning of each sequence in the nucleotide sequence files. For example, add the filename "SRR5935758_all.1" to the beginning of each sequence in the SRR5935758_all.1.bgc.fna file.
sed -i "s/>NODE/>SRR5935758_all.1_NODE/" SRR5935758_all.1.bgc.fna

#Merge all bgc.fna files with filenames added to the beginning of each sequence in the all_deepBGC_fna folder into a single file named all.bgc.fna.
find ./ -name "*.bgc.fna"| xargs cat > all.bgc.fna

# Use the SeqKit tool to extract the BGC gene set sequences according to the IDs in the first column of the deepBGC prediction results (all.bgc_2.tsv).
#Create a SeqKit virtual environment and install the SeqKit tool.
conda create --name seqkit
conda activate seqkit
conda install seqkit

#SeqKit requires the input fa file to be a compressed file, and the input ID file to be a txt file, so preprocess the input files first.
#Compress the all.bgc.fna into a gz file
gzip -c all.bgc.fna > all.bgc.fna.gz

#Keep only the first column of the all.bgc_2.tsv file as the ID file for the BGCs.
awk -F, '{print $1}' all.bgc_2.tsv > all.bgc_3.txt

#SeqKit extracts the BGC gene set sequences based on the IDs predicted by deepBGC.
zcat all.bgc.fna.gz | seqkit grep -f all.bgc_3.txt > picked_new.fa

#Compress the picked_new.fa into a gz file.
gzip -c picked_new.fa > picked_new.fa.gz

#Then split the BGC gene set sequences picked by BGCs ID (picked_new.fa.gz) into individual fa files.
seqkit split picked_new.fa.gz -i --id-regexp "^([\w]+)\-" -2

#Combine the BGC gene sets from the same sample to get the BGCs gene set for that sample.
cat *SRR5935758* > SRR5935758.fa
```
# Calculation of BGCs Relative Abundance

```bash
##Use salmon to calculate the relative abundance of BGCs.
#Create a salmon virtual environment and install the salmon tool.
conda create --name salmon
conda activate salmon
conda install salmon
#Build the index with salmon, where file_name.fa is the BGCs gene file from step 3.
salmon index -t file_name.fa -i transcripts_index
#Calculate the relative abundance of BGCs with salmon, the "quant.sf" file in the "file_name_transcripts_quant" folder contains the relative abundance of BGCs.
salmon quant -i transcripts_index  -l A -1 file_name_1.fastq -2 file_name_2.fastq -o file_name_transcripts_quant
#Convert the quant.sf format to txt format.
cat file_quant.sf > file_quant.txt
#Extract the first and fourth columns of the file and output to file_quant_1.txt. "file_quant_1.txt" is the relative abundance file for BGCs.
awk '{print $1,$4}' file_quant.txt > file_quant_1.txt
```
# Construction of BGCs Classification File

```bash
#The gut microbiome BGCs sequence (mag_humangut) can be replaced with other BGCs data you want to explore.
# Download the gut microbiome BGCs sequence (mag_humangut) from the BIG-FAM database and name it gut_BGC_backgroun.fasta, filter out sequences with aa less than 60 to become the background file (gut_BGC_backgroun_desolve.fasta).
awk 'BEGIN{OFS=FS="\t"}{if($0~/>/) name=$0 ;else seq[name]=seq[name]$0;}END{for(i in seq) {if(length(seq[i])>100) print i"\n"seq[i]}}' gut_BGC_backgroun.fasta > gut_BGC_backgroun_desolve.fasta
##Use the diamond tool to match the BGCs gene file and the background file to obtain the classification of BGCs.
#Build the diamond database.
diamond makedb --in gut_BGC_backgroun_desolve.fasta --db nr
#There are some sequences with trailing underscores in the BGCs gene file -- "file_name.fa" file, diamond will report an error and need to be deleted.
# Remove special characters (different numbers of underscores) from the nucleic acid file.
sed "/,/d" file_name.fa >  file_name_1.fa;sed 's/_____\+//g' file_name_1.fa >file_name_2.fa;sed 's/____\+//g' file_name_2.fa >file_name_3.fa;sed 's/___\+//g' file_name_3.fa >file_name_4.fa;sed 's/__\+//g' file_name_4.fa >file_name_5.fa
#Put the processed sequences into a new folder for subsequent processing.
find -name '*_5.fa' -exec cp {} /mnt/data1/ZhouJiaJia/HMP_IBD/HMP_1338/MG_1338/BGC_prediect_result_1338/deslove \;
#Diamond alignment.
diamond blastx --db nr -q file_name_5.fa -o file_name_fmt.txt 
# Filter the alignment results where the e-value is less than 10^-5 (the 11th column is the E value), and output to the file_name_fmt_1.txt file.
awk '{if($11<=1e-05) print $0}' file_name_fmt.txt > file_name_fmt_1.txt
#Extract the first two columns of the "file_name_fmt_1.txt" file to get the BGCs classification file.
awk '{print $1,$4}' file_name_fmt_1.txt > file_name_fmt_2.txt
```
# Enrichment of BGCs
#'file_quant_1.txt' is the relative abundance file of BGCs from step 4, and 'file_name_fmt_2.txt' is the classification file of BGCs from step 5.
gene_set<- read.csv('file1_name_fmt_2.txt',header = F,sep = " ")
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
ll = c()
ll2 = c()
for (i in 1:length(list)) {
  t = length(unlist(list[i]))
  if (t > 2) {
    ll = append(ll,i)
  }
  ll2 = append(ll2,t)
}
gsva_matrix<- gsva(as.matrix(uni_matrix), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
write.csv(gsva_matrix, "file2_BGC_score.csv")

uni_matrix<- read.table('file2_quant_1.txt',row.names=1,header = T,sep = " ")
gene_set<- read.csv('file2_name_fmt_2.txt',header = F,sep = " ")
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
ll = c()
ll2 = c()
for (i in 1:length(list)) {
  t = length(unlist(list[i]))
  if (t > 2) {
    ll = append(ll,i)
  }
  ll2 = append(ll2,t)
}
gsva_matrix<- gsva(as.matrix(uni_matrix), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
write.csv(gsva_matrix, "file2_BGC_score.csv")

uni_matrix<- read.table('file3_quant_1.txt',row.names=1,header = T,sep = " ")
gene_set<- read.csv('file3_name_fmt_2.txt',header = F,sep = " ")
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
ll = c()
ll2 = c()
for (i in 1:length(list)) {
  t = length(unlist(list[i]))
  if (t > 2) {
    ll = append(ll,i)
  }
  ll2 = append(ll2,t)
}
gsva_matrix<- gsva(as.matrix(uni_matrix), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
write.csv(gsva_matrix, "file3_BGC_score.csv")
#You can enter multiple loops.
# Write an R script to calculate the enrichment of BGCs.
Rscript mGSEA_enrichment.R