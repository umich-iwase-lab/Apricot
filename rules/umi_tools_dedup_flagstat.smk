#ALL.append(expand(OUTPUT_DIR + '/07-umi_tools_dedup_flagstat/{sample}.genome.umi_tools_dedup_flagstat.flagstat',
#                  sample=config['samples']))

rule umi_tools_dedup_flagstat:
    input:
        OUTPUT_DIR + '/07-umi_tools_dedup/{sample}.genome.umi_tools_dedup.bam'
    output:
        OUTPUT_DIR + '/07-umi_tools_dedup_flagstat/{sample}.genome.umi_tools_dedup_flagstat.flagstat',
    log:
        OUTPUT_DIR + '/07-umi_tools_dedup_flagstat/.log/{sample}.genome.umi_tools_dedup_flagstat.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/umi_tools_dedup_flagstat.{sample}.benchmark.txt'
    shell: '(samtools flagstat {input} > {output}) 2>&1 | tee {log}'
