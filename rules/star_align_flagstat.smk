ALL.append(expand(OUTPUT_DIR + '/04-star_align_flagstat/{sample}.star_align_flagstat.flagstat',
                  sample=config['samples']))

rule star_align_flagstat:
    input:
        OUTPUT_DIR + '/03-star_align/{sample}.star_align.bam'
    output:
        OUTPUT_DIR + '/04-star_align_flagstat/{sample}.star_align_flagstat.flagstat',
    log:
        OUTPUT_DIR + '/04-star_align_flagstat/.log/{sample}.star_align_flagstat.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/{sample}.star_align_flagstat.benchmark.txt'
    shell: '(samtools flagstat {input} > {output}) 2>&1 | tee {log}'
