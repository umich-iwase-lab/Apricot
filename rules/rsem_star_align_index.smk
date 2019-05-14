rule rsem_star_align_index:
    input:
        OUTPUT_DIR + '/rsem_star_align/{sample}.genome.bam'
    output:
        bam = OUTPUT_DIR + '/rsem_star_align/{sample}.genome.sorted.bam',
        bai = OUTPUT_DIR + '/rsem_star_align/{sample}.genome.sorted.bam.bai',
    threads: 12
    log:
        OUTPUT_DIR + '/rsem_star_align/.log/{sample}.rsem_star_align_index.log',
    benchmark:
        OUTPUT_DIR + '/benchmarks/rsem_star_align_index.{sample}.benchmark.txt'
    shell: '''(
    samtools sort -@ {threads} \
    -m 1G \
    -o {output.bam} \
    {input}

    samtools index {output.bam}) 2>&1 | tee {log}'''
