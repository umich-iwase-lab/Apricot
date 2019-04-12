rule samtools_sort_index_after_dedup:
    input:
        OUTPUT_DIR + '/06-umi_tools_dedup/{sample}.genome.umi_tools_dedup.bam'
    output:
        bam = OUTPUT_DIR + '/10-samtools_sort_index_after_dedup/{sample}.genome.umi_tools_dedup.sorted.bam',
        bai = OUTPUT_DIR + '/10-samtools_sort_index_after_dedup/{sample}.genome.umi_tools_dedup.sorted.bam.bai',
    threads: 12
    log:
        OUTPUT_DIR + '/10-sort_index_after_dedup/.log/{sample}.samtools_sort_index_after_dedup.log',
    benchmark:
        OUTPUT_DIR + '/benchmarks/samtools_sort_index_after_dedup.{sample}.benchmark.txt'
    shell: '''(
    samtools sort -@ {threads} \
    -m 1G \
    -o {output.bam} \
    {input}

    samtools index {output.bam}) 2>&1 | tee {log}'''
