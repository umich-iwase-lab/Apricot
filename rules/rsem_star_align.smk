_star_config = config['genome_reference']['star']

rule rsem_star_align:
    input:
        reads = lambda wildcards: expand(\
                OUTPUT_DIR + '/umi_tools_extract/processed.{sample}.{read}.fastq.gz',
                sample=wildcards.sample,
                read=READS,
                ),
        genomeParameters = join(\
                REFERENCE_DIR,
                _star_config['genome_dir'],
                'genomeParameters.txt'
                ),
    output:
        OUTPUT_DIR + '/rsem_star_align/{sample}.genes.results',
        OUTPUT_DIR + '/rsem_star_align/{sample}.isoforms.results',
        OUTPUT_DIR + '/rsem_star_align/{sample}.genome.bam',
        OUTPUT_DIR + '/rsem_star_align/{sample}.transcript.bam',
    log:
        OUTPUT_DIR + '/rsem_star_align/.log/{sample}.rsem_star_align.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/rsem_star_align.{sample}.benchmark.txt'
    threads: 12
    resources:
        mem_gb=30
    params:
        genomeDir = REFERENCE_DIR + '/' + _star_config['genome_dir'],
        outFileNamePrefix = OUTPUT_DIR + '/rsem_star_align/{sample}',
        paired_end = '--paired-end' if PAIRED_END else '',
    shell: '''(
STAR_PATH=$(dirname $(which STAR))
rsem-calculate-expression \
    {params.paired_end} \
    --star \
    --star-path $STAR_PATH \
    --star-gzipped-read-file \
    --star-output-genome-bam \
    -p {threads} \
    {input.reads} \
    {params.genomeDir}/RSEM_ref \
    {params.outFileNamePrefix}
mv {params.outFileNamePrefix}.STAR.genome.bam {params.outFileNamePrefix}.genome.bam
)2>&1 | tee {log}'''
