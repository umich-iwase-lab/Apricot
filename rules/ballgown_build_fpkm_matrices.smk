ALL.append([OUTPUT_DIR + '/10-ballgown_build_fpkm_matrices/fpkm-gene.txt',
            OUTPUT_DIR + '/10-ballgown_build_fpkm_matrices/fpkm-transcript.txt',
            ])

#TODO: please add logging on the R script using sink ala:
#https://github.com/snakemake-workflows/rna-seq-star-deseq2/blob/master/scripts/deseq2.R#L1
rule ballgown_build_fpkm_matrices:
    input:
        ballgown_dir = expand(OUTPUT_DIR + '/08-stringtie/ballgown/{sample}',
                              sample=config['samples']),
    output:
        fpkm_gene = OUTPUT_DIR + '/10-ballgown_build_fpkm_matrices/fpkm-gene.txt',
        fpkm_transcript = OUTPUT_DIR + '/10-ballgown_build_fpkm_matrices/fpkm-transcript.txt',
    benchmark:
        OUTPUT_DIR + '/benchmarks/ballgown_build_fpkm_matrices.benchmark.txt'
    params:
        output_dir = OUTPUT_DIR + '/10-ballgown_build_fpkm_matrices',
        input_dir= OUTPUT_DIR + '/08-stringtie/ballgown',
    script: '../scripts/ballgown_build_fpkm_matrices.R'
