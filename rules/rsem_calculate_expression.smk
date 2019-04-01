rule rsem_calculate_expression:
    input:
        OUTPUT_DIR + '/05-umi_tools_dedup/{sample}.umi_tools_dedup.bam'
    output:
        OUTPUT_DIR + '/08-rsem_calculate_expression/{sample}.genes.results',
        OUTPUT_DIR + '/08-rsem_calculate_expression/{sample}.isoforms.results',
    log:
        OUTPUT_DIR + '/08-rsem_calculate_expression/.log/{sample}.rsem_calculate_expression.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/{sample}.rsem_calculate_expression.txt'
    params:
        genomeDir = REFERENCE_DIR + '/' + _star_config['genome_dir'],
        outFileNamePrefix = OUTPUT_DIR + '/08-rsem_calculate_expression/{sample}',
    shell:'''(
rsem-calculate-expression --alignments \
--paired-end \
{input} \
{params.genomeDir}/RSEM_ref \
{params.outFileNamePrefix}
    )2>&1 | tee {log}'''
