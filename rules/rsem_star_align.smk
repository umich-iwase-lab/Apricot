_star_config = config['genome_reference']['star']

#ALL.append(expand(OUTPUT_DIR + '/03-rsem_star_align/{sample}.transcript.sorted.bam.bai',
#               sample=config['samples']))
rule rsem_star_align:
    input:
        reads = lambda wildcards: expand(\
                OUTPUT_DIR + '/01-umi_tools_extract/processed.{sample}.{read}.fastq.gz',
                sample=wildcards.sample,
                read=['R1', 'R2'],
                ),
        genomeParameters = join(\
                REFERENCE_DIR,
                _star_config['genome_dir'],
                'genomeParameters.txt'
                ),
    output:
        OUTPUT_DIR + '/03-rsem_star_align/{sample}.genes.results',
        OUTPUT_DIR + '/03-rsem_star_align/{sample}.isoforms.results',
        OUTPUT_DIR + '/03-rsem_star_align/{sample}.STAR.genome.bam',
        OUTPUT_DIR + '/03-rsem_star_align/{sample}.transcript.bam',
        #OUTPUT_DIR + '/03-rsem_star_align/{sample}.transcript.sorted.bam',
        #OUTPUT_DIR + '/03-rsem_star_align/{sample}.transcript.sorted.bam.bai',
    log:
        OUTPUT_DIR + '/03-rsem_star_align/.log/{sample}.rsem_star_align.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/rsem_star_align.{sample}.benchmark.txt'
    threads: 12
    resources:
        mem_gb=30
    params:
        genomeDir = REFERENCE_DIR + '/' + _star_config['genome_dir'],
#        sjdbGTFfile = REFERENCE_DIR + '/' + config['genome_reference']['gtf'],
        outFileNamePrefix = OUTPUT_DIR + '/03-rsem_star_align/{sample}',
        paired_end = '--paired-end' if config['sequencing_parameters']['paired'] else '',
        #star_bam_file = '{sample}.Aligned.sortedByCoord.out.bam',
        #output_bam_file = '{sample}.star_align.bam'
    shell: '''(
rsem-calculate-expression \
{params.paired_end} \
--star \
--star-path ~/miniconda3/envs/apricot/bin \
--star-gzipped-read-file \
--star-output-genome-bam \
-p {threads} \
{input.reads} \
{params.genomeDir}/RSEM_ref \
{params.outFileNamePrefix}
    )2>&1 | tee {log}
    '''
#--sort-bam-by-coordinate \
#     shell: '''(
# STAR \
#     --twopassMode Basic \
#     --genomeDir {params.genomeDir} \
#     --runThreadN {threads} \
#     --readFilesIn {input.reads} \
#     --readFilesCommand 'gunzip -c' \
#     --outFileNamePrefix {params.outFileNamePrefix} \
#     --sjdbGTFfile {params.sjdbGTFfile} \
#     --runMode alignReads \
#     --outSAMtype BAM SortedByCoordinate \
#     --limitBAMsortRAM 64000000000 \
#     --outFilterType BySJout \
#     --outFilterMultimapNmax 20 \
#     --alignSJoverhangMin 8 \
#     --alignSJDBoverhangMin 1 \
#     --alignIntronMin 20 \
#     --outFilterMismatchNmax 999 \
#     --outFilterScoreMinOverLread 0 \
#     --outFilterMatchNminOverLread 0 \
#     --alignMatesGapMax 1000000 \
#     --alignIntronMax 1000000 \
#     --chimSegmentMin 20
#
# pushd {OUTPUT_DIR}/03-star_align/
# ln --symbolic {params.star_bam_file} {params.output_bam_file}
# popd
# ) 2>&1 | tee {log}'''
