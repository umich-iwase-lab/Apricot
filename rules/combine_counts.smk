ALL.append(expand(OUTPUT_DIR + '/04-combine_counts/{feature}_{metric}.txt',
                  feature=['gene','isoform'],
                  metric=['expected_count', 'FPKM', 'TPM']))
            # OUTPUT_DIR + '/04-combine_counts/gene_fpkm.txt',
            # OUTPUT_DIR + '/04-combine_counts/gene_tpm.txt',
            # OUTPUT_DIR + '/04-combine_counts/isoform_count.txt',
            # OUTPUT_DIR + '/04-combine_counts/isoform_fpkm.txt',
            # OUTPUT_DIR + '/04-combine_counts/isoform_tpm.txt',
            # ])
rule combine_counts:
    input:
        genes = expand(OUTPUT_DIR + '/03-rsem_star_align/{sample}.genes.results',
                      sample=config['samples']),
        isoforms = expand(OUTPUT_DIR + '/03-rsem_star_align/{sample}.isoforms.results',
                      sample=config['samples']),
    output:
        gene = OUTPUT_DIR + '/04-combine_counts/gene_{metric}.txt',
        isoform = OUTPUT_DIR + '/04-combine_counts/isoform_{metric}.txt',

        # gene_count = OUTPUT_DIR + '/04-combine_counts/gene_count.txt',
        # gene_fpkm = OUTPUT_DIR + '/04-combine_counts/gene_fpkm.txt',
        # gene_tpm = OUTPUT_DIR + '/04-combine_counts/gene_tpm.txt'
        # isoform_count = OUTPUT_DIR + '/04-combine_counts/isoform_count.txt',
        # isoform_fpkm = OUTPUT_DIR + '/04-combine_counts/isoform_fpkm.txt',
        # isoform_tpm = OUTPUT_DIR + '/04-combine_counts/isoform_tpm.txt',

    benchmark:
        OUTPUT_DIR + '/benchmarks/combine_counts_{metric}.benchmark.txt'
    log:
        OUTPUT_DIR + '/04-combine_counts/.log/combine_counts_{metric}.log'
    params:
        metric = '{metric}',
        gene_count_script = srcdir('../scripts/combine.py'),
    shell:'''(
python {params.gene_count_script} --input_path 'outputs/03-rsem_star_align/*.genes.results' \
--output_file {output.gene} \
-c {params.metric} \
--id_columns gene_id

python {params.gene_count_script} --input_path 'outputs/03-rsem_star_align/*.isoforms.results' \
--output_file {output.isoform} \
-c {params.metric} \
--id_columns transcript_id,gene_id


)2>&1 | tee {log}
    '''
