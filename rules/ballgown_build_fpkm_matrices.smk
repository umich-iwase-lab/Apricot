ALL.append(expand(OUTPUT_DIR + '/02-fastqc_seq/processed.{sample}.{read}_fastqc.html',
                  sample=config['samples'],
                  read=['R1', 'R2']))
rule ballgown_build_fpkm_matrices:
    input:
        #ballgon dirs
    output:
        #gene and transcript fpkm matrices
    benchmark:
        'benchmarks/ballgown_build_fpkm_matrices.benchmark.txt'
    params:
        output_dir = ?
    shell:
        '''
        Rscript scripts/ballgown_build_fpkm_matrices.R \
            --input_dir {input.?} \
            --output_dir {params.output_dir}
        '''
