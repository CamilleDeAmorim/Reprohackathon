Samples=["SRX7080612","SRX7080613","SRX7080614","SRX7080615","SRX7080616","SRX7080617"]

rule all: 
    input:
        expand("results/01_raw_data/{SRA_id}.fastq.gz", SRA_id=Samples),
        expand("results/02_Trimming_results/{SRA_id}_trimmed.fq.gz", SRA_id=Samples),
        ("results/03_Reference_Genome/reference.fasta"), 
        ("results/04_Genome_Annotation/reference.gff")
# Rule that download the FASTQ files required for the analysis.
rule fasterq_dump:
    output:
        "results/01_raw_data/{SRA_id}.fastq.gz"
    container:
        "images/fasterq-dump.img"
    log:
        "logs/fasterq-dump/{SRA_id}.log"
    threads: 60
    shell:
        """
        fasterq-dump --threads {threads} --progress {wildcards.SRA_id} -O results/01_raw_data/ &> {log}
        gzip -f results/01_raw_data/{wildcards.SRA_id}.fastq
        """
# Rule that performs the trimming of FASTQ files.
rule trim_galore:
    input:
        "results/01_raw_data/{SRA_id}.fastq.gz"
    output: 
        trimmed_fastq="results/02_Trimming_results/{SRA_id}_trimmed.fq.gz",
        reports="results/02_Trimming_results/{SRA_id}_trimming_report.txt"
    container:
        "images/TrimGalor.img"
    log:
        "logs/TrimGalor/{SRA_id}_trimmed.log"
    threads: 60
    shell:
        """
        trim_galore -q 20 --phred33 --length 25 {input} -o results/02_Trimming_results/ &> {log}
        mv results/02_Trimming_results/{wildcards.SRA_id}.fastq.gz_trimming_report.txt {output.reports}

        """
# Rule that downloads the reference genome in FASTA format.
rule reference_genome:
    output:
        "results/03_Reference_Genome/reference.fasta"
    shell:
        """
        wget -q -O  results/03_Reference_Genome/reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
        """
# Rule that downloads the genome annotation.
rule genome_annotation :
    output:
        "results/04_Genome_Annotation/reference.gff"
    log: 
        "logs/Genome_Annotation"
    shell:
        """
        wget -O {output} "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"  &> {log}       
        """

