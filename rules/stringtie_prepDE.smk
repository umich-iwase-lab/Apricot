ALL.append([OUTPUT_DIR + '/09-stringtie_prepDE/gene_count_matrix.stringtie_prepDE.csv',
            OUTPUT_DIR + '/09-stringtie_prepDE/transcript_count_matrix.stringtie_prepDE.csv',
            ])

rule stringtie_prepDE:
    input:
        gtfs = expand(OUTPUT_DIR + '/08-stringtie/{sample}.gtf',
                      sample=config['samples']),
        prep_input = expand(OUTPUT_DIR + '/08-stringtie/{sample}_prepDE_input.txt',
                            sample=config['samples']),
    output:
        gene_counts = OUTPUT_DIR + '/09-stringtie_prepDE/gene_count_matrix.stringtie_prepDE.csv',
        transcript_counts = OUTPUT_DIR + '/09-stringtie_prepDE/transcript_count_matrix.stringtie_prepDE.csv',
    log:
        OUTPUT_DIR + '/09-stringtie_prepDE/.log/stringtie_prepDE.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/stringtie_prepDE.benchmark.txt'
    params:
        length = config['sequencing_parameters']['read_length'],
        prep_input = ' '.join(expand(OUTPUT_DIR + '/08-stringtie/{sample}_prepDE_input.txt',
                                     sample=config['samples'])),
        prep_config = OUTPUT_DIR + '/09-stringtie_prepDE/prepDE_config.txt',
    shell:
        '''(
cat {params.prep_input} > {params.prep_config}
prepDE.py \
    -i {params.prep_config} \
    -g {output.gene_counts} \
    -t {output.transcript_counts} \
    -l {params.length}
) 2>&1 | tee {log}'''
