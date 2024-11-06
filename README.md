# mGSEA: metagenomic single sample Gene Set Enrichment Analysis  

mGSEA is a novel tool developed based on the principle of ssGSEA.
mGSEA is a tool that can be applied to the analysis of metagenomic or metatranscriptomic data.
mGSEA can calculate the enrichment scores of BGCs in individual metagenomic samples.
For mGSEA analysis, two key files are required: the abundance data of BGCs and the classification information of BGCs.

## 第一步_conda安装deepBGC

```bash
#安装deepBGC工具（用于预测肠道菌群中的BGC。若有更好的工具，可替换该工具），deepBGC需要在Conda创建的环境中使用。
#添加bioconda源和conda-forge。
conda config --add channels bioconda
conda config --add channels conda-forge
#创建deepBGC虚拟环境，安装依赖工具。
conda create -n deepbgc python=3.7 hmmer prodigal
#激活deepBGC虚拟环境。
conda activate deepbgc
#安装deepBGC。
conda install deepbgc
#下载训练好的模型和Pfam数据库。
deepbgc download
#查看下载的依赖和模型。
deepbgc info
```
#第二步_deepBGC工具预测BGC及预测结果的筛选

```bash
#deepBGC工具预测BGC。“file_scaffolds.fasta为组装后的contig文件”，（加“--prodigal-meta-mode”可以检测短contigs中的基因）。
deepbgc pipeline --prodigal-meta-mode  file_scaffolds.fasta
#每个样本，deepBGC预测后都会生成一个xx.bgc.tsv后缀的结果文件，该文件包含多条序列的预测结果，为避免批量处理时，混淆结果，先在每个xx.bgc.tsv文件的第一列添加上样本名称。"^"表示第一列 "sed -i"表示直接在文件中修改，不打印出来。
for file in xx.bgc.tsv; do sed -i "s/^/$file/" "$file"; done
#合并所有bgc.tsv后缀的deepBGC预测结果文件，方便批量处理，生成一个名为all.bgc.tsv的汇总文件。
find ./ -name "*.bgc.tsv"| xargs cat > all.bgc.tsv
#删掉all.bgc.tsv文件第18列（product_class列）为空值的行，输出新结果到all.bgc_1.tsv中。
awk -F, '$11 != ""' all.bgc.tsv > all.bgc_1.tsv
#筛选all.bgc_1.tsv文件第12列（deepbgc_score列）阈值大于0.5的结果，输出到all.bgc_2.tsv文件中。all.bgc_2.tsv文件即为deepBGC预测的BGC结果，包含BGC基因集的名称及可能的预测分类。
awk -F, '$3 > 0.5' all.bgc_1.tsv > all.bgc_2.tsv
```
#第三步_deepBGC基因集的构建

```bash
#将deepBGC预测结果中的所有的genbank文件汇总到一个名为“all_deepBGC_gbk”的新文件夹，便于批量处理。
find -name '*.gbk' -exec cp {} /mnt/data1/ZhouJiaJia/HMP_IBD/HMP_1338/MG_1338/all_deepBGC_gbk \;
#查看all_deepBGC_gbk文件中有哪些genbank文件。
cd /mnt/data1/ZhouJiaJia/HMP_IBD/HMP_1338/MG_1338/all_deepBGC_gbk && ls
#用pick_fa.py脚本从genbank文件中提取核苷酸序列。
python pick_fa.py
#此处列举了部分例子。
```
#python
#撰写名为pick_fa.py的python脚本从genbank文件中提取核酸序列，"SRR5935758_all.1.bgc.gbk"、"SRR5935758_all.2.bgc.gbk"、"SRR5946523_all.1.bgc.gbk".....genbank文件。
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
#将从genbank文件中提取核苷酸序列文件汇总到一个名为“all_deepBGC_fna”的新文件夹，便于批量处理。
find -name '*.bgc.fna' -exec cp {} /mnt/data1/ZhouJiaJia/HMP_IBD/HMP_1338/MG_1338/all_deepBGC_fna \;
#查看all_deepBGC_fna文件中有哪些核苷酸序列文件。
cd /mnt/data1/ZhouJiaJia/HMP_IBD/HMP_1338/MG_1338/all_deepBGC_fna && ls
#为每个核苷酸序列文件中的每条序列前加上文件名。例如在SRR5935758_all.1.bgc.fna文件中的每条序列前加上文件名“SRR5935758_all.1”。
sed -i "s/>NODE/>SRR5935758_all.1_NODE/" SRR5935758_all.1.bgc.fna
#将all_deepBGC_fna文件夹中，每条序列前加上文件名后的所有bgc.fna后缀的文件合并为all.bgc.fna。
find ./ -name "*.bgc.fna"| xargs cat > all.bgc.fna
#用SeqKit工具按第一步中筛选后的deepBGC预测结果（all.bgc_2.tsv）中第一列的ID提取对应的BGC基因集序列。
#创建SeqKit虚拟环境，安装SeqKit工具。
conda create --name seqkit
conda activate seqkit
conda install seqkit
#seqkit工具要求输入的fa文件需为压缩文件，输入的ID文件为txt文件，需要先对待输入文件进行处理。
#将all.bgc.fna压缩为gz文件
gzip -c all.bgc.fna > all.bgc.fna.gz
#只保留all.bgc_2.tsv文件的第一列作为BGCs的ID文件
awk -F, '{print $1}' all.bgc_2.tsv > all.bgc_3.txt
#SeqKit根据deepBGC预测的BGC的ID提取BGC的基因集序列。
zcat all.bgc.fna.gz | seqkit grep -f all.bgc_3.txt > picked_new.fa
#将picked_new.fa压缩为gz文件。
gzip -c picked_new.fa > picked_new.fa.gz
#再按id将按BGCs ID挑选后的BGCs基因集序列（picked_new.fa.gz）拆分为单个fa文件。
seqkit split picked_new.fa.gz -i --id-regexp "^([\w]+)\-" -2
#把同一样本的BGC基因集合并到一起，即为该样本的BGCs基因集。
cat *SRR5935758* > SRR5935758.fa
```
#第四步_BGCs相对丰度的计算

```bash
##salmon计算BGCs的相对丰度。
#创建salmon虚拟环境，安装salmon工具。
conda create --name salmon
conda activate salmon
conda install salmon
#salmon构建索引,file_name.fa是第三步中的BGCs基因文件。
salmon index -t file_name.fa -i transcripts_index
#salmon计算BGCs的相对丰度,"file_name_transcripts_quant"文件夹中的quant.sf就是包含BGCs相对丰度的文件。
salmon quant -i transcripts_index  -l A -1 file_name_1.fastq -2 file_name_2.fastq -o file_name_transcripts_quant
#将quant.sf格式转为txt格式。
cat file_quant.sf > file_quant.txt
#提取文件的第一列和第四列，输出到file_quant_1.txt。“file_quant_1.txt”即为BGCs的相对丰度文件。
awk '{print $1,$4}' file_quant.txt > file_quant_1.txt
```
# 第五步_构建BGCs分类文件

```bash
#肠道菌群BGCs序列（mag_humangut）可替换为其他自己想探究的BGCs数据。
#从BIG-FAM数据库下载肠道菌群BGCs序列（mag_humangut）命名为gut_BGC_backgroun.fasta，过滤掉其中aa小于60的序列，变成背景库文件（gut_BGC_backgroun_desolve.fasta）。
awk 'BEGIN{OFS=FS="\t"}{if($0~/>/) name=$0 ;else seq[name]=seq[name]$0;}END{for(i in seq) {if(length(seq[i])>100) print i"\n"seq[i]}}' gut_BGC_backgroun.fasta > gut_BGC_backgroun_desolve.fasta
##用diamond工具比对BGCs的基因文件和背景库文件以获取BGCs的分类。
#diamond建库。
diamond makedb --in gut_BGC_backgroun_desolve.fasta --db nr
BGCs的基因文件——“file_name.fa”文件中，有些序列末尾会有一些下划线，diamond会报错，需要删掉。
#删掉核酸文件中的特殊字符(不同数量的下划线)。
sed "/,/d" file_name.fa >  file_name_1.fa;sed 's/_____\+//g' file_name_1.fa >file_name_2.fa;sed 's/____\+//g' file_name_2.fa >file_name_3.fa;sed 's/___\+//g' file_name_3.fa >file_name_4.fa;sed 's/__\+//g' file_name_4.fa >file_name_5.fa
#把处理好的序列放入新的文件夹，以便后续处理。
find -name '*_5.fa' -exec cp {} /mnt/data1/ZhouJiaJia/HMP_IBD/HMP_1338/MG_1338/BGC_prediect_result_1338/deslove \;
#diamond比对。
diamond blastx --db nr -q file_name_5.fa -o file_name_fmt.txt 
#筛选比对结果中，e值小于10-5次方的结果（第11列是E value值），输出到file_name_fmt_1.txt文件中。
awk '{if($11<=1e-05) print $0}' file_name_fmt.txt > file_name_fmt_1.txt
#提取“file_name_fmt_1.txt”文件的前2列，即为BGCs的分类文件。
awk '{print $1,$4}' file_name_fmt_1.txt > file_name_fmt_2.txt
```
#第六步_BGCs的富集
#R语言
#'file_quant_1.txt'是第四步中BGCs的相对丰度文件，'file_name_fmt_2.txt'是第五步中BGCs的分类文件。
uni_matrix<- read.table('file1_quant_1.txt',row.names=1,header = T,sep = " ")
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
#可以输入多个循环。
```bash
#编写R脚本，计算BGCs的富集程度。
Rscript mGSEA_enrichment.R
```
