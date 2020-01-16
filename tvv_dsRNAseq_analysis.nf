#!/usr/bin/env nextflow

// Setup default parameters if not provided any
params.sampleName = "123414"
params.reads = "input_reads/*sub*{1,2}.fastq.gz"
params.outputDirectory = "output/"
params.phiX = "genomes/phiX.fasta"
params.tvvDir = "genomes/"

// Read in the parameters as usable variables
sample_name = params.sampleName
output_directory = params.outputDirectory
tvv_directory = params.tvvDir

// Read in the sample's reads files
paired_reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )

// Read in phiX genome for depletion
phiX_ch = Channel.fromPath( params.phiX, checkIfExists: true )

// Read in the TVV reference genomes for mapping & binning
tvv_ch = Channel.fromPath( "${tvv_directory}/tvv*.fasta", checkIfExists: true )

// Split it into one input channel per set of sequences
tvv_ch.into { tvv1_ch; tvv2_ch; tvv3_ch; tvv4_ch; tvv5_ch; tvv_sats_ch }

// Trim adapters from the reads
process trimAdapters {
    publishDir "${output_directory}/analysis/01_adapter_trimmed/", mode: "copy"

    input:
    tuple val(sampleID), file(reads) from paired_reads_ch

    output:
    file "*val*fq.gz" into trimmed_reads

    """
    trim_galore \
        --paired \
        --stringency 5 \
        --quality 20 \
        --fastqc \
        --gzip \
        --basename $sample_name \
        $reads
    """
}

// Remove any reads that map to phiX (spiked in during RNA-seq library prep for sequencing
process depletePhiX {
    publishDir "${output_directory}/analysis/02_phiX_depleted/", mode: "copy"

    input:
    file trimmed_reads
    file phiX_genome from phiX_ch

    output:
    file "*phiX_depleted*.R1.fq.gz" into phiX_depleted_forward_reads
    file "*phiX_depleted*.R2.fq.gz" into phiX_depleted_reverse_reads
    file "*phiX.stats.txt"
    file "*phiX.counts.txt"
    """
    # Build BWA index out of the reference
    bwa index \
    -p phiX_index \
    $phiX_genome

    # Perform the mapping
    bwa mem \
    phiX_index \
    $trimmed_reads > phiX_mapped_sam

    # Get summary stats of the mapping
    samtools flagstat phiX_mapped_sam > ${sample_name}.phiX.stats.txt

    # Retrieve phiX-depleted reads & sort output
    samtools view -f 4 -bh phiX_mapped_sam | samtools sort > ${sample_name}.phiX_depleted.sorted.bam

    # Save how many phiX reads were in the sequencing data
    samtools idxstats ${sample_name}.phiX_depleted.sorted.bam > \
        ${sample_name}.phiX.counts.txt

    # Print that phiX count to the screen
    echo "name   length   mapped_reads   unmapped_reads"
    cat ${sample_name}.phiX.counts.txt

    # Convert the phiX-depleted BAM back to individual FASTQ files
    samtools fastq \
    -1 ${sample_name}.phiX_depleted.R1.fq \
    -2 ${sample_name}.phiX_depleted.R2.fq \
    -s ${sample_name}.phiX_depleted.singletons.fq \
    ${sample_name}.phiX_depleted.sorted.bam

    # Compress the phiX-depleted fastqs
    gzip ${sample_name}.phiX_depleted.*fq
    """
}

process mapToTVV {
    publishDir "${output_directory}/analysis/03_binned_reads/fastq/", mode: "copy"

    input:
    file forward_reads from phiX_depleted_forward_reads
    file reverse_reads from phiX_depleted_reverse_reads
    file tvv1 name "tvv1.fasta" from tvv1_ch
    file tvv2 name "tvv2.fasta" from tvv2_ch
    file tvv3 name "tvv3.fasta" from tvv3_ch
    file tvv4 name "tvv4.fasta" from tvv4_ch
    file tvv5 name "tvv5.fasta" from tvv5_ch
    file tvv_satellites name "tvv-dsRNA-satellites.fasta" from tvv_sats_ch

    output:
    file "${sample_name}_tvv*.fq" into binned_fastq_ch

    """
    bbsplit.sh \
        in1=$forward_reads \
        in2=$reverse_reads \
        ref="${tvv1},${tvv2},${tvv3},${tvv4},${tvv5},${tvv_satellites}" \
        basename="${sample_name}_%.fq" \
        outu=unmapped_reads \
        ambig2=best
    """
}

process fastqToFasta {
    publishDir "${output_directory}/analysis/03_binned_reads/fasta/", mode: "copy"

    input:
    file fastq from binned_fastq_ch.flatten()

    output:
    file "*fasta"

    """
    seqtk seq -A $fastq > ${fastq.baseName}.fasta
    """
}
