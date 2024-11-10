# apptainer installation
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer-suid


apptainer build images/fasterq-dump.img container/Singularity.fasterq-dump
apptainer build images/TrimGalor.img container/Singularity.TrimGalor
apptainer build images/cutadapt_v1.11.img container/Singularity.cutadapt_v1.11
apptainer build images/bowtie_v0.12.7_samtools.img container/Singularity.bowtie_v0.12.7_samtools
apptainer build images/featureCounts_v1.4.6-p3.img container/Singularity.featureCounts_v1.4.6-p3


