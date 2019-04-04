ALL.append(expand(OUTPUT_DIR + '/09-combine_counts_after_dedup/{feature}_{metric}.txt',
                  feature=['gene','isoform'],
                  metric=['expected_count', 'FPKM', 'TPM']))

rule combine_counts_after_dedup:
    input:
        genes = expand(OUTPUT_DIR + '/08-rsem_calculate_expression_using_dedup/{sample}.genes.results',
                      sample=config['samples']),
        isoforms = expand(OUTPUT_DIR + '/08-rsem_calculate_expression_using_dedup/{sample}.isoforms.results',
                      sample=config['samples']),
    output:
        gene = OUTPUT_DIR + '/09-combine_counts_after_dedup/gene_{metric}.txt',
        isoform = OUTPUT_DIR + '/09-combine_counts_after_dedup/isoform_{metric}.txt',
    benchmark:
        OUTPUT_DIR + '/benchmarks/combine_counts_after_dedup_{metric}.benchmark.txt'
    log:
        OUTPUT_DIR + '/09-combine_counts_after_dedup/.log/combine_counts_after_dedup_{metric}.log'
    params:
        metric = '{metric}',
        gene_count_script = srcdir('../scripts/combine.py'),
    shell:'''(
python {params.gene_count_script} --input_path 'outputs/08-rsem_calculate_expression_using_dedup/*.genes.results' \
--output_file {output.gene} \
-c {params.metric} \
--id_columns gene_id

python {params.gene_count_script} --input_path 'outputs/08-rsem_calculate_expression_using_dedup/*.isoforms.results' \
--output_file {output.isoform} \
-c {params.metric} \
--id_columns transcript_id,gene_id

)2>&1 | tee {log}
    '''
