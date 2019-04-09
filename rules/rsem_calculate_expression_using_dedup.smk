# ALL.append(OUTPUT_DIR + '/08-rsem_calculate_expression_using_dedup/{sample}.genes.results')
# ALL.append(OUTPUT_DIR + '/08-rsem_calculate_expression_using_dedup/{sample}.isoforms.results')
#ALL.append(expand(OUTPUT_DIR + '/08-rsem_calculate_expression_using_dedup/{sample}.{feature}.results',
#                  sample=config['samples'], feature=['genes','isoforms']))
rule rsem_calculate_expression_using_dedup:
    input: #fix to include genome
        OUTPUT_DIR + '/07-convert_dedup_for_rsem/{sample}.transcript.convert_dedup_for_rsem.bam',
    output:
        OUTPUT_DIR + '/08-rsem_calculate_expression_using_dedup/{sample}.genes.results',
        OUTPUT_DIR + '/08-rsem_calculate_expression_using_dedup/{sample}.isoforms.results',
    log:
        OUTPUT_DIR + '/08-rsem_calculate_expression_using_dedup/.log/{sample}.rsem_calculate_expression_using_dedup.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/{sample}.rsem_calculate_expression_using_dedup.txt'
    params:
        genomeDir = REFERENCE_DIR + '/' + _star_config['genome_dir'],
        outFileNamePrefix = OUTPUT_DIR + '/08-rsem_calculate_expression_using_dedup/{sample}',
    shell:'''(
rsem-calculate-expression --alignments \
--paired-end \
--no-bam-output \
{input} \
{params.genomeDir}/RSEM_ref \
{params.outFileNamePrefix}
    )2>&1 | tee {log}'''
