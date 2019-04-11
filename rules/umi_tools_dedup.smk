rule umi_tools_dedup:
    input:
        # bam = OUTPUT_DIR + '/03-rsem_star_align/{sample}.transcript.sorted.bam',
        # bai = OUTPUT_DIR + '/03-rsem_star_align/{sample}.transcript.sorted.bam.bai',
        bam = OUTPUT_DIR + '/04A-exclude_non_primary/{sample}.{type}.exclude_non_primary.bam',
        bai = OUTPUT_DIR + '/04A-exclude_non_primary/{sample}.{type}.exclude_non_primary.bam.bai',
        #bam = OUTPUT_DIR + '/03-rsem_star_align/{sample}.rsem.STAR.genome.sorted.bam',
        #bai = OUTPUT_DIR + '/03-rsem_star_align/{sample}.rsem.STAR.genome.sorted.bam.bai',
    output:
        OUTPUT_DIR + '/06-umi_tools_dedup/{sample}.{type}.umi_tools_dedup.bam',
    log:
        OUTPUT_DIR + '/06-umi_tools_dedup/.log/{sample}.{type}.umi_tools_dedup.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/umi_tools_dedup.{type}.{sample}.benchmark.txt'
    params:
        paired_end = '--paired' if config['sequencing_parameters']['paired'] else '',
    shell: '''(
umi_tools dedup \
    -I {input.bam} \
    {params.paired_end} \
    -S {output}\
) 2>&1 | tee {log}'''
