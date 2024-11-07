Samples = ["SRX7080612", "SRX7080613", "SRX7080614", "SRX7080615", "SRX7080616", "SRX7080617"]

# Rule that downloads all necessary files.
rule all: 
    input:
        expand("results/01_raw_data/{SRA_id}.fastq.gz", SRA_id=Samples),
        expand("results/02_Trimming_results/{SRA_id}_trimmed.fq.gz", SRA_id=Samples),
        "results/03_Reference_Genome/reference.fasta", 
        "results/04_Genome_Annotation/reference.gff"

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
    threads: 40
    shell:
        """
        trim_galore -q 20 --phred33 --length 25 {input} -o results/02_Trimming_results/ &> {log}
        
        """
#mv results/02_Trimming_results/{wildcards.SRA_id}.fastq.gz_trimming_report.txt {output.reports}
# Rule that downloads the reference genome in FASTA format.
rule reference_genome:
    output:
        "results/03_Reference_Genome/reference.fasta"
    shell:
        """
        wget -q -O {output} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"
        """

# Rule that downloads the genome annotation.
rule genome_annotation:
    output:
        "results/04_Genome_Annotation/reference.gff"
    log: 
        "logs/Genome_Annotation"
    shell:
        """
        wget -O results/04_Genome_Annotation/reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"  &> {log}       
        """

# Rule that creates the genome index.
rule genome_index:
    input:
        "results/03_Reference_Genome/reference.fasta"
    output:
        "results/03_Reference_Genome/reference.1.ebwt",
        "results/03_Reference_Genome/reference.2.ebwt",
        "results/03_Reference_Genome/reference.3.ebwt",
        "results/03_Reference_Genome/reference.4.ebwt",
        "results/03_Reference_Genome/reference.rev.1.ebwt",
        "results/03_Reference_Genome/reference.rev.2.ebwt"
    shell: 
        """
        bowtie-build {input} results/03_Reference_Genome/reference
        """
