rule fastqc_align:
    input:
        OUTPUT_DIR + '/07-umi_tools_dedup/{sample}.genome.umi_tools_dedup.bam',
    output:
        OUTPUT_DIR + '/10-fastqc_align/{sample}_fastqc.html',
    log:
        OUTPUT_DIR + '/10-fastqc_align/.log/{sample}.fastqc_align.log'
    benchmark:
        OUTPUT_DIR + '/benchmarks/fastqc_align.{sample}.benchmark.txt'
    params:
        symlink = OUTPUT_DIR + '/10-fastqc_align/{sample}.bam',
        fastqc_dir= OUTPUT_DIR+ '/10-fastqc_align/',
    shell: '''(
python -c 'import os, os.path; os.symlink(os.path.relpath("{input}", os.path.dirname("{params.symlink}")), "{params.symlink}")'
fastqc {params.symlink} -o {params.fastqc_dir}
) 2>&1 | tee {log}'''
