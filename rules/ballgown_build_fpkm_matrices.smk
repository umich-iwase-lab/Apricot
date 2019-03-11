ALL.append([OUTPUT_DIR + '/10-ballgown_build_fpkm_matrices/gene_fpkms.txt',
            OUTPUT_DIR + '/10-ballgown_build_fpkm_matrices/iso_fpkms.txt',
            ])

rule ballgown_build_fpkm_matrices:
    input:
        ballgown_dir = expand(OUTPUT_DIR + '/08-stringtie/ballgown/{sample}',
                              sample=config['samples']),
    output:
        fpkm_gene = OUTPUT_DIR + '/10-ballgown_build_fpkm_matrices/iso_fpkms.txt',
        fpkm_transcript = OUTPUT_DIR + '/10-ballgown_build_fpkm_matrices/gene_fpkms.txt',
    benchmark:
        OUTPUT_DIR + '/benchmarks/ballgown_build_fpkm_matrices.benchmark.txt'
    params:
        output_dir = OUTPUT_DIR + '/10-ballgown_build_fpkm_matrices',
        input_dir= OUTPUT_DIR + '/08-stringtie/ballgown',
    script: '../scripts/ballgown_build_fpkm_matrices.R'
