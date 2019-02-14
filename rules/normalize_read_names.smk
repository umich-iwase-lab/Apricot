rule normalize_read_names:
    input:
        reads = lambda wildcards: expand(INPUT_DIR + '/{read_name}',
                                         read_name=config['samples'][wildcards.sample]),
    output:
        reads = expand(OUTPUT_DIR + '/00-normalize_read_names/{{sample}}.{read}.fastq.gz',
                       read=['R1', 'R2']),
    run:
        from os import symlink
        from os.path import dirname
        from os.path import relpath
        for in_read, out_read in zip(input.reads, output.reads):
            symlink(relpath(in_read, dirname(out_read)), out_read)
