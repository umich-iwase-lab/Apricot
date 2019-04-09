ALL.append(expand(OUTPUT_DIR + '/04-rsem_star_align_flagstat/{sample}.rsem_star_align_flagstat.flagstat',
                  sample=config['samples']))

rule rsem_star_align_flagstat:
    input:
        OUTPUT_DIR + '/03-rsem_star_align/{sample}.STAR.genome.bam'
    output:
        OUTPUT_DIR + '/04-rsem_star_align_flagstat/{sample}.rsem_star_align_flagstat.flagstat',
    log:
        OUTPUT_DIR + '/04-rsem_star_align_flagstat/.log/{sample}.rsem_star_align_flagstat.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/rsem_star_align_flagstat.{sample}.benchmark.txt'
    shell: '(samtools flagstat {input} > {output}) 2>&1 | tee {log}'
