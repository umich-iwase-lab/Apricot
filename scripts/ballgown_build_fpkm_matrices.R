#TODO Add doc
input_dir = snakemake@params[["input_dir"]]
output_dir = snakemake@params[["output_dir"]]

sample_paths = list.files(input_dir, full.names = TRUE)

if(length(sample_paths) == 0) {
    stop(sprintf('No sample directories listed in %s', input_dir))
}

library(yaml)
library(stringr)
library(ballgown)

gene_fpkm_file = sprintf('%s/fpkm-gene.txt', output_dir)
iso_fpkm_file = sprintf('%s/fpkm-transcript.txt', output_dir)

#######################################
message('Reading ballgown input')
bg_data = ballgown(samples = sample_paths, meas = 'all')

###################
# Pull out gene-level FPKMs
gene_fpkms = gexpr(bg_data)
# Retain the id column (gene symbols)
gene_fpkms = cbind(id = rownames(gene_fpkms), data.frame(gene_fpkms))

###################
# Pull out transcript-level FPKMs
iso_fpkms = texpr(bg_data, meas = 'all')

# Create a dictionary between t_id (internal to ballgown), gene_id (symbol),
# the t_name (transcript RefSeq name), and transcript location / strand
# NOTE: This information is determined from the gtf files input to ballgown()
# which are in turn determined based on the gtf given in the config file.
gtf_dict = iso_fpkms[, c('t_id', 'gene_id', 't_name', 'chr', 'start', 'end', 'strand')]

# Retain identifying information for the transcripts
iso_col_keep = c('t_id', 'gene_id', 't_name', 'chr', 'start', 'end', 'strand', grep('FPKM', colnames(iso_fpkms), value = T))
iso_fpkms = iso_fpkms[, iso_col_keep]

message(sprintf('Writing gene FPKMs to %s', gene_fpkm_file))
write.table(gene_fpkms, file = gene_fpkm_file, sep = '\t', row.names = F, col.names = T, quote = F)

message(sprintf('Writing isoform FPKMs to %s', iso_fpkm_file))
write.table(iso_fpkms, file = iso_fpkm_file, sep = '\t', row.names = F, col.names = T, quote = F)
