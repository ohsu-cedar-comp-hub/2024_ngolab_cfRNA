library(HGNChelper)

#current versions contain corrected gene names using this script.
#you likely do not need to run this script. It is being included for posterity.
#set this to true if you do want to run this
fix_gene_names <- FALSE
genes_tocheck <- c("1-Dec", "1-Mar", "1-Sep", "10-Mar", "10-Sep", 
                        "11-Mar", "11-Sep", "12-Sep", "14-Sep", "2-Mar", 
                        "2-Sep", "3-Mar", "3-Sep", "4-Mar", "4-Sep", 
                        "5-Mar", "5-Sep", "6-Mar", 
                        "6-Sep", "7-Mar", "7-Sep",  "8-Mar", 
                        "8-Sep", "9-Mar", "9-Sep")
corrected_name <- HGNChelper::checkGeneSymbols(genes_tocheck)
rownames(corrected_name) <- corrected_name$x

#now we correct the gene names that were messed up for file one
genecount_file <- "./Analysis_Scripts/pdac_genecount.csv.gz"
pdac_genecount<- read.csv(genecount_file,
                          row.names=1)
corrected_rowname<- rownames(pdac_genecount)
corrected_rowname <- ifelse(corrected_rowname %in% genes_tocheck,
                            corrected_name[corrected_rowname, "Suggested.Symbol"],
                            corrected_rowname) %>% str_replace(".\\/\\/\\/.", "_or_due_to_excel_")
rownames(pdac_genecount) <- corrected_rowname

if (fix_gene_names){
  write.csv(pdac_genecount, gzfile(genecount_file), row.names = T)
  print("written")
}

#now correcting for the normed file as well
normed_count_file <- "./Analysis_Scripts/normalized_pdac_genecount.csv.gz"
pdac_genecount_normed <- read.csv(normed_count_file,
                          row.names=1)
corrected_rowname<- rownames(pdac_genecount_normed)
corrected_rowname <- ifelse(corrected_rowname %in% genes_tocheck,
                            corrected_name[corrected_rowname, "Suggested.Symbol"],
                            corrected_rowname) %>% str_replace(".\\/\\/\\/.", "_or_due_to_excel_")
rownames(pdac_genecount_normed) <- corrected_rowname

if (fix_gene_names){
  write.csv(pdac_genecount_normed, gzfile(normed_count_file), row.names = T)
  print("written")
}



