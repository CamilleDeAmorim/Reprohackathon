if [ ! -d workflow_execution ] ; then mkdir workflow_execution; fi
if [ ! -d workflow_execution/data ] ; then 
    mkdir workflow_execution/data; 
    apptainer run container/fasterq-dump.img --threads 12 --progress PRJNA586837 --outdir workflow_execution/data/ ;
fi
if [ ! -d workflow_execution/trimmed ] ; then mkdir workflow_execution/trimmed ; fi
apptainer run container/trimgalore.img -q 20 --phred33 --length 25 workflow_execution/data/PRJNA586837.fastq --output_dir workflow_execution/trimmed/

wget -q -O workflow_execution/data/reference.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta"

apptainer exec container/bowtie.img bowtie-build workflow_execution/data/reference.fasta workflow_execution/data/genome_index

if [ ! -d workflow_execution/mapping ] ; then mkdir workflow_execution/mapping ; fi
apptainer run container/bowtie.img -p 8 -S workflow_execution/data/genome_index workflow_execution/trimmed/PRJNA586837_trimmed.fq | apptainer run container/samtools.img sort -@ 8 > workflow_execution/mapping/PRJNA586837.bam
apptainer run container/samtools.img index workflow_execution/mapping/PRJNA586837.bam

wget -O workflow_execution/data/reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"

mkdir -p workflow_execution/featureCounts
apptainer run container/featurecounts.img --extraAttributes Name -t gene -g ID -F GTF -T 8 -a workflow_execution/data/reference.gff -o workflow_execution/featureCounts/counts.txt workflow_execution/mapping/PRJNA586837.bam