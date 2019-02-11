ALL.append(expand(OUTPUT_DIR + '/06-umi_tools_dedup_flagstat/{sample}.umi_tools_dedup_flagstat.flagstat',
                  sample=config['samples']))

rule umi_tools_dedup_flagstat:
    input:
        OUTPUT_DIR + '/05-umi_tools_dedup/{sample}.umi_tools_dedup.bam'
    output:
        OUTPUT_DIR + '/06-umi_tools_dedup_flagstat/{sample}.umi_tools_dedup_flagstat.flagstat',
    log:
        OUTPUT_DIR + '/06-umi_tools_dedup_flagstat/.log/{sample}.umi_tools_dedup_flagstat.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/{sample}.umi_tools_dedup_flagstat.benchmark.txt'
    shell: '(samtools flagstat {input} > {output}) 2>&1 | tee {log}'
