ALL.append(OUTPUT_DIR + '/90-multiqc/multiqc_report.html')

rule multiqc:
    input: expand(OUTPUT_DIR + '/02-fastqc_seq/processed.{sample}.{read}_fastqc.html', sample=config['samples'], read=READS),
           expand(OUTPUT_DIR + '/04-rsem_star_align_flagstat/{sample}.rsem_star_align_flagstat.flagstat',sample=config['samples']),
           expand(OUTPUT_DIR + '/08-umi_tools_dedup_flagstat/{sample}.genome.umi_tools_dedup_flagstat.flagstat',sample=config['samples']),
           expand(OUTPUT_DIR + '/10-fastqc_align/{sample}_fastqc.html',sample=config['samples']),
    output: OUTPUT_DIR + '/90-multiqc/multiqc_report.html'
    log: OUTPUT_DIR + '/90-multiqc/.log/multiqc_report.log'
    benchmark: OUTPUT_DIR + '/benchmarks/multiqc.benchmark.txt'
    params:
        multiqc_dir= OUTPUT_DIR + '/90-multiqc/',
        config_file = srcdir('../config/multiqc_config.yaml'),
    shell: '''
(multiqc \
    {OUTPUT_DIR} \
    -o {params.multiqc_dir} \
    --config {params.config_file}
) 2>&1 | tee {log}'''
