library(data.table)

# outputs

# inputs
circ_dir <- snakemake@input[["circ_dir"]]

# params
contrast <- snakemake@params[["contrast"]]
print(contrast)
first <- contrast[1]
second <- contrast[2]
print(first)
group <- snakemake@params[["linear_model"]]
sample_id <- snakemake@params[["sample_id"]]

metadata <- snakemake@params[["samples"]]
md <- read.delim(file=metadata, stringsAsFactors = FALSE)

#dir.create(paste0("ciri2_output/", first,"-vs-",second), showWarnings = FALSE)

first_list <- md[which(md[,group] == first), sample_id]
second_list <- md[which(md[,group] == second), sample_id]
print(first_list)
namelist <- append(first_list,second_list)
newlist <- list()
for (i in namelist) {
namelist[which(namelist==i)] <- paste0("ciri2_output/",first,"-vs-", second,"/ciri2out_bwa_map_", i,".txt")
}

firstlen <- length(first_list)
secondlen <- length(second_list)
lengths <- append(firstlen, secondlen)
lengths <- paste0(lengths, collapse=",")
#namelist <- paste0(newlist, collapse=",")
write(namelist, file = paste0("ciri2_output/", first, "-vs-", second,"/", first, "-vs-", second,"_namelist.txt"))
write(lengths, file = paste0("ciri2_output/", first, "-vs-", second,"/", first, "-vs-", second,"_namelengths.txt"))
