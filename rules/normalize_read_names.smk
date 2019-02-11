rule normalize_read_names:
    input:
        read1 = lambda wildcards: INPUT_DIR + '/' + config['samples'][wildcards.sample][0],
        read2 = lambda wildcards: INPUT_DIR + '/' + config['samples'][wildcards.sample][1],
    output:
        read1 = OUTPUT_DIR + '/00-normalize_read_names/{sample}.R1.fastq.gz',
        read2 = OUTPUT_DIR + '/00-normalize_read_names/{sample}.R2.fastq.gz',
    shell: '''
ln -fsr {input.read1} {output.read1}
ln -fsr {input.read2} {output.read2}'''
