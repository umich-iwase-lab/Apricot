ALL.extend(expand(OUTPUT_DIR + '/combine_counts_after_dedup/{feature}_{metric}_dedup.txt',
                  feature=['gene','isoform'],
                  metric=['expected_count', 'FPKM', 'TPM']))

_INPUT_DIR = OUTPUT_DIR + '/rsem_calculate_expression_using_dedup/'

rule combine_counts_after_dedup:
    input:
        genes = expand(_INPUT_DIR + '{sample}.genes.results',
                       sample=config['samples']),
        isoforms = expand(_INPUT_DIR + '{sample}.isoforms.results',
                          sample=config['samples']),
    output:
        gene = OUTPUT_DIR + '/combine_counts_after_dedup/gene_{metric}_dedup.txt',
        isoform = OUTPUT_DIR + '/combine_counts_after_dedup/isoform_{metric}_dedup.txt',
    benchmark:
        OUTPUT_DIR + '/benchmarks/combine_counts_after_dedup_{metric}.benchmark.txt'
    log:
        OUTPUT_DIR + '/combine_counts_after_dedup/.log/combine_counts_after_dedup_{metric}.log'
    params:
        metric = '{metric}',
        gene_count_script = srcdir('../scripts/combine.py'),
    shell:'''(
python {params.gene_count_script} \
    --input_path '{_INPUT_DIR}*.genes.results' \
    --output_file {output.gene} \
    -c {params.metric} \
    --id_columns gene_id
python {params.gene_count_script} \
    --input_path '{_INPUT_DIR}*.isoforms.results' \
    --output_file {output.isoform} \
    -c {params.metric} \
    --id_columns transcript_id,gene_id
)2>&1 | tee {log}
    '''
