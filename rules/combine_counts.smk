ALL.extend(expand(OUTPUT_DIR + '/combine_counts/{feature}_{metric}.txt',
                  feature=['gene','isoform'],
                  metric=['expected_count', 'FPKM', 'TPM']))

_INPUT_DIR = OUTPUT_DIR + '/rsem_star_align/'

rule combine_counts:
    input:
        genes = expand(_INPUT_DIR + '{sample}.genes.results',
                      sample=config['samples']),
        isoforms = expand(_INPUT_DIR +'{sample}.isoforms.results',
                      sample=config['samples']),
    output:
        gene = OUTPUT_DIR + '/combine_counts/gene_{metric}.txt',
        isoform = OUTPUT_DIR + '/combine_counts/isoform_{metric}.txt',
    benchmark:
        OUTPUT_DIR + '/benchmarks/combine_counts_{metric}.benchmark.txt'
    log:
        OUTPUT_DIR + '/combine_counts/.log/combine_counts_{metric}.log'
    params:
        metric = '{metric}',
        gene_count_script = srcdir('../scripts/combine.py'),
    shell:'''(
python {params.gene_count_script} --input_path '{_INPUT_DIR}*.genes.results' \
--output_file {output.gene} \
-c {params.metric} \
--id_columns gene_id

python {params.gene_count_script} --input_path '{_INPUT_DIR}/*.isoforms.results' \
--output_file {output.isoform} \
-c {params.metric} \
--id_columns transcript_id,gene_id


)2>&1 | tee {log}
    '''
