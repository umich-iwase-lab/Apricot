rule umi_tools_dedup:
    input:
        bam = OUTPUT_DIR + '/samtools_exclude_non_primary/{sample}.{type}.samtools_exclude_non_primary.bam',
        bai = OUTPUT_DIR + '/samtools_exclude_non_primary/{sample}.{type}.samtools_exclude_non_primary.bam.bai',
    output:
        OUTPUT_DIR + '/umi_tools_dedup/{sample}.{type}.umi_tools_dedup.bam',
    log:
        OUTPUT_DIR + '/umi_tools_dedup/.log/{sample}.{type}.umi_tools_dedup.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/umi_tools_dedup.{type}.{sample}.benchmark.txt'
    params:
        paired_end = '--paired' if PAIRED_END else '',
    shell: '''(
umi_tools dedup \
    -I {input.bam} \
    {params.paired_end} \
    -S {output}\
) 2>&1 | tee {log}'''
