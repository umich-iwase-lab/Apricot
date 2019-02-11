rule star_align:
    input:
        OUTPUT_DIR + '/01-umi_tools_extract/processed.{sample}.R1.fastq.gz',
        OUTPUT_DIR + '/01-umi_tools_extract/processed.{sample}.R2.fastq.gz',
    output:
        OUTPUT_DIR + '/03-star_align/{sample}.star_align.bam'
    log:
        OUTPUT_DIR + '/03-star_align/.log/{sample}.star_align.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/{sample}.star_align.benchmark.txt'
    threads: 12
    params:
        genomeDir = REFERENCE_DIR + '/' + config['genome_reference']['star_genome_dir'],
        sjdbGTFfile = REFERENCE_DIR + '/' + config['genome_reference']['gtf'],
        outFileNamePrefix = OUTPUT_DIR + '/03-star_align/{sample}.',
        star_bam_file = '{sample}.Aligned.sortedByCoord.out.bam',
        output_bam_file = '{sample}.star_align.bam'
    shell: '''(
STAR \
--genomeDir {params.genomeDir} \
--runThreadN {threads} \
--readFilesIn {input} \
--readFilesCommand 'gunzip -c' \
--outFileNamePrefix {params.outFileNamePrefix} \
--sjdbGTFfile {params.sjdbGTFfile} \
--runMode alignReads \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 64000000000 \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--alignIntronMin 20 \
--outFilterMismatchNmax 999 \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0 \
--alignMatesGapMax 1000000 \
--alignIntronMax 1000000 \
--chimSegmentMin 20

pushd {OUTPUT_DIR}/03-star_align/
ln --symbolic {params.star_bam_file} {params.output_bam_file}
popd
) 2>&1 | tee {log}'''
