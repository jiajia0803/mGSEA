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