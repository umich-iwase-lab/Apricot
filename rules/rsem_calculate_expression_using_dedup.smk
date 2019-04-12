rule rsem_calculate_expression_using_dedup:
    input:
        OUTPUT_DIR + '/09-convert_dedup_for_rsem/{sample}.transcript.convert_dedup_for_rsem.bam',
    output:
        OUTPUT_DIR + '/11-rsem_calculate_expression_using_dedup/{sample}.genes.results',
        OUTPUT_DIR + '/11-rsem_calculate_expression_using_dedup/{sample}.isoforms.results',
    log:
        OUTPUT_DIR + '/11-rsem_calculate_expression_using_dedup/.log/{sample}.rsem_calculate_expression_using_dedup.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/{sample}.rsem_calculate_expression_using_dedup.txt'
    params:
        genomeDir = REFERENCE_DIR + '/' + _star_config['genome_dir'],
        outFileNamePrefix = OUTPUT_DIR + '/11-rsem_calculate_expression_using_dedup/{sample}',
        paired_end = '--paired-end' if config['sequencing_parameters']['paired'] else '',
    shell:'''(
rsem-calculate-expression --alignments \
    {params.paired_end} \
    --no-bam-output \
    {input} \
    {params.genomeDir}/RSEM_ref \
    {params.outFileNamePrefix}
)2>&1 | tee {log}'''
