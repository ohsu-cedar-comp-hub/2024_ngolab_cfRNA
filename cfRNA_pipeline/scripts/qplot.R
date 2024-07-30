library("data.table")
library("qvalue")

stats_table <- snakemake@input[['stats_table']]
cat(sprintf(c('stats table: ', stats_table, '\n')))

list_stats_table <- snakemake@input[['list_stats_table']]
cat(sprintf(c('stats table: ', list_stats_table, '\n')))

ERCC_stats_table <- snakemake@input[['ERCC_stats_table']]
cat(sprintf(c('stats table: ', ERCC_stats_table, '\n')))

qplot <- snakemake@output[['qplot']]
cat(sprintf(c('Qvalue Output: ', qplot, '\n')))

qhist <- snakemake@output[['qhist']]
cat(sprintf(c('Qvalue hist Output: ', qhist, '\n')))

out_table = snakemake@output[['table']]
cat(sprintf(c('Summary results table', out_table,'\n')))

list_qplot <- snakemake@output[['list_qplot']]
cat(sprintf(c('Qvalue Output: ', list_qplot, '\n')))

list_qhist <- snakemake@output[['list_qhist']]
cat(sprintf(c('Qvalue hist Output: ', list_qhist, '\n')))

list_out_table = snakemake@output[['list_table']]
cat(sprintf(c('Summary results table', list_out_table,'\n')))

ERCC_qplot <- snakemake@output[['ERCC_qplot']]
cat(sprintf(c('Qvalue Output: ', ERCC_qplot, '\n')))

ERCC_qhist <- snakemake@output[['ERCC_qhist']]
cat(sprintf(c('Qvalue hist Output: ', ERCC_qhist, '\n')))

ERCC_out_table = snakemake@output[['ERCC_table']]
cat(sprintf(c('Summary results table', ERCC_out_table,'\n')))

stats_frame = read.table(stats_table, row.names=1, sep='\t', check.names=F)

qobj = qvalue(p=stats_frame$pvalue, fdr.level=T)

stats_frame$qvalues = qobj$qvalues
stats_frame$lfdr = qobj$lfdr
write.table(as.data.frame(stats_frame), file=out_table, quote=FALSE, sep='\t')

pdf(qplot)
plot(qobj)
dev.off()

pdf(qhist)
hist(qobj)
dev.off()

############ List normalized

stats_frame = read.table(list_stats_table, row.names=1, sep='\t', check.names=F)

qobj = qvalue(p=stats_frame$pvalue, fdr.level=T)

stats_frame$qvalues = qobj$qvalues
stats_frame$lfdr = qobj$lfdr
write.table(as.data.frame(stats_frame), file=out_table, quote=FALSE, sep='\t')

pdf(list_qplot)
plot(qobj)
dev.off()

pdf(list_qhist)
hist(qobj)
dev.off()

############### ERCC Normalized

stats_frame = read.table(ERCC_stats_table, row.names=1, sep='\t', check.names=F)

qobj = qvalue(p=stats_frame$pvalue, fdr.level=T)

stats_frame$qvalues = qobj$qvalues
stats_frame$lfdr = qobj$lfdr
write.table(as.data.frame(stats_frame), file=out_table, quote=FALSE, sep='\t')

pdf(ERCC_qplot)
plot(qobj)
dev.off()

pdf(ERCC_qhist)
hist(qobj)
dev.off()




