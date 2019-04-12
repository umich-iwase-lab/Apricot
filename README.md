# Apricot
Iwase lab RNA-seq pipeline

Apricot is a bulk RNA-Seq pipeline that accepts in-line barcoded, paried fastq files and produces sample alignments, count matrices, sample bigwigs, and an array of QC reports. It uses Snakemake as an automation framework and has been tested in several Unix/Linux environments. It uses STAR (wrapped by RSEM) for alignment and as such will requires peak memory on the order of 30Gb RAM to align using mm9/mm10.


Quickstart
----------
### Setup new machine or new user
  - Clone Apricot from github
    ```
    cd ~; git clone <github url>
    ```
  - If necessary, install Miniconda following instructions here:
    https://conda.io/projects/conda/en/latest/user-guide/install/index.html
  - Build and activate conda environment
    ```
    conda env create -f Apricot/envs/apricot.yml #first time only
    conda activate apricot
    ```
  - Review versions
    ```
    # brief list of key program versions
    cat Apricot/envs/apricot.yml
    # full list of versions (including dependencies)
    conda list
    ```
### Add a new genome reference
  ```
  cd reference
  mkdir GRCm38_vM21; cd GRCm38_vM21
  ```
  - Download genome fasta and gtf files
    https://www.gencodegenes.org/mouse/
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/GRCm38.primary_assembly.genome.fa.gz
  - Add ERCC spike-ins if desired
    - Download "ERCC Controls Annotation: ERCC RNA Spike-In Control Mixes" from Thermo Fisher
      https://www.thermofisher.com/order/catalog/product/4456739
    - Append sequences to fasta file
      ```
      wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095047.txt
      # gunzip original fasta/gtf as necessary
      cat GRCm38_vM21.fa ercc.fa | dos2unix > GRCm38_vM21_ercc.fa

    - Append gtf entries to gtf file
      ```
      cat cms_095047.txt \
      | awk 'NR>1 {printf("%s\tercc\texon\t1\t%i\t.\t+\t.\tgene_id\"G%s\";gene_version \"1\";gene_name \"%s\";gene_source \ \"ercc\";gene_biotype \"ercc\";transcript_id \"%s\";\n", \
      $1, length($5)-1, $1, $1, $1)}' > ercc.gtf

      cat GRCm38_vM21.gtf ercc.gtf | dos2unix > GRCm38_vM21_ercc.gtf

- Setup project
  ```
  mkdir -p project_A/inputs
  cd project_A
  ```
  - copy sample fastq files into project_A/inputs

- Setup configfile
  - Copy config/config.template.yaml to project_A/config.project_A.yaml
  - Edit to review/revise genome, directories, sequencing parameters
  - Add sample names and fastq names

- Run workflow
  - Setup screen and activate the conda environment
    ```screen -S project_A
    conda activate apricot
    ```
  - Show rule graph
    ```
    snakemake \
        --snakefile ~/Apricot/Snakefile \
        --configfile config.project_A.yaml \
        --rulegraph | dot -Tpdf > rulegraph.pdf
    ```
  - Dry-run to check rules and jobs that will be executed:
    ```
    snakemake \
        --snakefile ~/Apricot/Snakefile \
        --configfile config.project_A.yaml \
        --dryrun --quiet
    ```

  - Run, setting cores and memory
    ```
    snakemake \
        --snakefile ../Apricot/Snakefile \
        --configfile config.project_vallianatos.mm10_ercc.yaml \
        --resources mem_gb=105 \
        --cores 36 \
        -p
    ```

- Key outputs
  - Main outputs
    - combine_counts: the FPKM, TPM, count matrices prior to deduplication
    - combine_counts_after_dedup: FPKM, TPM, count matrices useful as input to
      DESeq2
    - deeptools_bigwig: bigwig files for each sample
    - multiqc/multiqc_report.html: summary of pipeline results
  - Logs
    - Snakemake run logs: all details as printed to screen. Note if running jobs
      in parallel, job details will sometimes be interleaved.
    ```
    cd {your project dir}
    cd .snakemake/log
    ls #list of logs ordered by run timestamp
    less `ls -tr | tail -1` #log details of last run
    ```
    - Rule logs: each rule logs details of its execution to a .log directory
      alongside main outputs
    - Benchmarks: all jobs log timing and memory use to the benchmarks
      directory alongside main outputs
      ```
      awk 'BEGIN {OFS="\t"} NR==1 {print "filename", $0} FNR>1 {print FILENAME, $0}' benchmarks/*.txt
      ```
      - Slowest 20 jobs:
        ```
        awk 'BEGIN {OFS="\t"} NR==1 {print "filename", $0} FNR>1 {print FILENAME, $0}' benchmarks/*.txt | sort -k2,2nr | head -20
        ```

### Tool Citations

- fastqc: Andrews S. (2010). FastQC: a quality control tool for high throughput
  sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
- multiqc: MultiQC: Summarize analysis results for multiple tools and samples
  in a single report. Philip Ewels, Måns Magnusson, Sverker Lundin and Max Käller. Bioinformatics (2016). doi: 10.1093/bioinformatics/btw354. PMID: 27312411
- samtools: Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N.,
  et al. (2009). The sequence alignment / map format and SAMtools. Bioinformatics 25 2078–2079.
- snakemake: Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable
  bioinformatics workflow engine”. Bioinformatics 2012.
- RSEM: Accurate transcript quantification from RNA-Seq data with or without a
  reference genome. Li B, Dewey C. BMC Bioinformatics, 2011, Volume 12.
  doi: 10.1186/1471-2105-12-323
- deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis.  
  Ramírez, F, et al.  Nucleic Acids Research (2016). doi:10.1093/nar/gkw257.
- STAR: Dobin, A et al. STAR: ultrafast universal RNA-seq aligner.
  Bioinformatics29, 15–21 (2013).
- UMITools: Smith TS, Heger A, Sudbery I (2017) UMI-tools: modelling sequencing
  errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Res 27: 491 – 499.
- The Sequence alignment/map (SAM) format and SAMtools. Li, H. et al. and 1000
  Genome Project Data Processing Subgroup (2009) Bioinformatics, 25, 2078-9. [PMID: 19505943]

### Basic orientation

- Snakemake
  - https://snakemake.readthedocs.io/en/stable/index.html
  - https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html

- Conda
  - https://conda.io/docs/
