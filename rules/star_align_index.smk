rule star_align_index:
    input:
        OUTPUT_DIR + '/03-star_align/{sample}.star_align.bam'
    output:
        OUTPUT_DIR + '/03-star_align/{sample}.star_align.bam.bai',
    log:
        OUTPUT_DIR + '/03-star_align/.log/{sample}.star_align_index.log',
    benchmark:
        OUTPUT_DIR + '/benchmarks/star_align_index.{sample}.benchmark.txt'
    shell: '(samtools index {input}) 2>&1 | tee {log}'
