ALL.append(expand(OUTPUT_DIR + '/11-deeptools/{sample}.bw',sample=config['samples']))

rule deeptools:
    input: OUTPUT_DIR + '/03-star_align/{sample}.star_align.bam'
    output: OUTPUT_DIR + '/11-deeptools/{sample}.bw'
    benchmark: OUTPUT_DIR + '/benchmarks/deeptools.{sample}.benchmark.txt'
    log:
        OUTPUT_DIR + '/11-deeptools/.log/{sample}.deeptools.log'
    shell: '(bamCoverage -b {input} -o {output}) 2>&1 | tee {log}'
