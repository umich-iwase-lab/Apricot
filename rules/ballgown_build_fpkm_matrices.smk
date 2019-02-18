ALL.append([OUTPUT_DIR + '/10-ballgown/fpkm-gene.txt',
            OUTPUT_DIR + '/10-ballgown/fpkm-transcript.txt',
            ])

rule ballgown_build_fpkm_matrices:
    input:
        ballgown_dir = expand(OUTPUT_DIR + '/08-stringtie/ballgown/{sample}',
                              sample=config['samples']),
    output:
        fpkm_gene = OUTPUT_DIR + '/10-ballgown/fpkm-gene.txt',
        fpkm_transcript = OUTPUT_DIR + '/10-ballgown/fpkm-transcript.txt',
    benchmark:
        OUTPUT_DIR + '/benchmarks/ballgown_build_fpkm_matrices.benchmark.txt'
    params:
        output_dir = 'foo'
    shell:
        '''
        touch {output.fpkm_gene}
        touch {output.fpkm_transcript}
#        Rscript scripts/ballgown_build_fpkm_matrices.R \
#            --input_dir {input.ballgown_dir} \
#            --output_dir {params.output_dir}
        '''
