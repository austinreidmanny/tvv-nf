#!/usr/bin/env nextflow

// Setup default parameters if not provided any
params.sample = "sample"
params.reads = "input/123414*one-percent*R{1,2}.fastq"
params.outputDirectory = "output/" + params.sample
params.phix = "genomes/phiX.fasta"
params.tvvDirectory = "genomes/"
params.diamondDB = "/n/data1/hms/mbib/nibert/austin/diamond/ncbi-viruses_tvv5.dmnd"
params.blockSize = "2"
params.threads = "4"

// Read in the sample's reads files
paired_reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )

// Read in phiX genome for depletion
phiX_ch = Channel.fromPath( params.phix, checkIfExists: true )

// Read in the TVV reference genomes for mapping & binning
tvv1_ch = Channel.fromPath(params.tvvDir + "/tvv1.fasta", checkIfExists: true)
            .map { file -> tuple(file.getSimpleName(), file) }
tvv2_ch = Channel.fromPath(params.tvvDir + "/tvv2.fasta", checkIfExists: true)
            .map { file -> tuple(file.getSimpleName(), file) }
tvv3_ch = Channel.fromPath(params.tvvDir + "/tvv3.fasta", checkIfExists: true)
            .map { file -> tuple(file.getSimpleName(), file) }
tvv4_ch = Channel.fromPath(params.tvvDir + "/tvv4.fasta", checkIfExists: true)
            .map { file -> tuple(file.getSimpleName(), file) }
tvv5_ch = Channel.fromPath(params.tvvDir + "/tvv5.fasta", checkIfExists: true)
            .map { file -> tuple(file.getSimpleName(), file) }
tvv_sats_ch = Channel.fromPath(params.tvvDir + "/tvv-dsRNA-satellites.fasta", checkIfExists: true)
            .map { file -> tuple(file.getSimpleName(), file) }

// Trim adapters from the reads
process trimAdapters {
    publishDir "${params.outputDirectory}/analysis/01_adapter_trimmed/", mode: "copy"

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
        --basename $params.sample \
        $reads
    """
}

// Remove any reads that map to phiX (spiked in during RNA-seq library prep for sequencing
process depletePhiX {
    publishDir "${params.outputDirectory}/analysis/02_phiX_depleted/", mode: "copy"

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
    samtools flagstat phiX_mapped_sam > ${params.sample}.phiX.stats.txt

    # Retrieve phiX-depleted reads & sort output
    samtools view -f 4 -bh phiX_mapped_sam | samtools sort > ${params.sample}.phiX_depleted.sorted.bam

    # Save how many phiX reads were in the sequencing data
    samtools idxstats ${params.sample}.phiX_depleted.sorted.bam > \
        ${params.sample}.phiX.counts.txt

    # Print that phiX count to the screen
    echo "name   length   mapped_reads   unmapped_reads"
    cat ${params.sample}.phiX.counts.txt

    # Convert the phiX-depleted BAM back to individual FASTQ files
    samtools fastq \
    -1 ${params.sample}.phiX_depleted.R1.fq \
    -2 ${params.sample}.phiX_depleted.R2.fq \
    -s ${params.sample}.phiX_depleted.singletons.fq \
    ${params.sample}.phiX_depleted.sorted.bam

    # Compress the phiX-depleted fastqs
    gzip ${params.sample}.phiX_depleted.*fq
    """
}

process mapToTVV {
    publishDir "${params.outputDirectory}/analysis/03_binned_reads/fastq/", mode: "copy"

    input:
    file forward_reads from phiX_depleted_forward_reads
    file reverse_reads from phiX_depleted_reverse_reads
    tuple tvv_species, file(tvv1) from tvv1_ch
    tuple tvv_species, file(tvv2) from tvv2_ch
    tuple tvv_species, file(tvv3) from tvv3_ch
    tuple tvv_species, file(tvv4) from tvv4_ch
    tuple tvv_species, file(tvv5) from tvv5_ch
    tuple tvv_species, file(tvv_satellites) from tvv_sats_ch

    output:
    file "${params.sample}_*fq" into binned_tvv, binned_tvv_fastq

    """
    bbsplit.sh \
        in1=$forward_reads \
        in2=$reverse_reads \
        ref="${tvv1},${tvv2},${tvv3},${tvv4},${tvv5},${tvv_satellites}" \
        basename="${params.sample}_%.fq" \
        outu=unmapped_reads \
        ambig2=best
    """
}

process fastqToFasta {
    publishDir "${params.outputDirectory}/analysis/03_binned_reads/fasta/", mode: "copy"

    input:
    file fastq from binned_tvv_fastq.flatten()

    output:
    file "${params.sample}_*.fasta"

    """
    seqtk seq -A $fastq > "${fastq.simpleName}.fasta"
    """
}


process deNovoAssembly {

    publishDir path: "${params.outputDirectory}/analysis/04_denovoassembly",
               pattern: "transcripts.fasta",
               mode: "copy",
               saveAs: { filename -> "${reads.getSimpleName()}.${filename}" }

    input:
    file reads from binned_tvv.flatten()

    output:
    tuple file("transcripts.fasta"), file(reads) into tvv_contigs

    // Build contigs with rnaSPAdes; if it cannot build any contigs for that TVV species, output an empty file
    """
    # Assemble viral reads into viral genomes or partial contigs
    rnaspades.py -s $reads -o ./ || touch transcripts.fasta
    """
}

// Map the TVV-reads to the denovo-assembled-contigs with BWA
process mapReadsToContigs {

    // Save the BAM file with the name of the TVV species
    publishDir path: "${params.outputDirectory}/analysis/05_refinement/mapping/",
               pattern: "reads_mapped_to_contigs.sorted.bam",
               mode: "copy",
               saveAs: { filename -> "${reads.getSimpleName()}.${filename}" }

    // Save the mapping-statistics file with the name of the TVV species
    publishDir path: "${params.outputDirectory}/analysis/05_refinement/mapping/",
              pattern: "reads_mapped_to_contigs.sorted.stats",
              mode: "copy",
              saveAs: { filename -> "${reads.getSimpleName()}.${filename}" }

    // Only read in the files for TVV species where rnaSPAdes could actually construct contigs
    input:
    tuple file(contigs), file(reads) from tvv_contigs.filter { it.get(0).size() > 0 }

    output:
    tuple file("reads_mapped_to_contigs.sorted.bam"), file(contigs), file(reads) into mapped_bam
    file "reads_mapped_to_contigs.stats"

    """
    # Create a BWA index of the contigs
    bwa index \
    -p contigs_index \
    $contigs

    # Map the reads to the contigs
     bwa mem \
     contigs_index \
     $reads > reads_mapped_to_contigs.sam

    # Get summary stats of the mapping
    samtools flagstat \
    reads_mapped_to_contigs.sam > \
    reads_mapped_to_contigs.stats

    # Remove unmapped reads and sort output (save only contig-mapped reads
    samtools view -F 4 -bh reads_mapped_to_contigs.sam | \
    samtools sort - > reads_mapped_to_contigs.sorted.bam

    # Remove (very large) uncompressed sam file
    rm reads_mapped_to_contigs.sam

    echo 'complete!'
    """
}


// With the reads mapped to their denovo-assemblies, refine the assembly by
// finding any mismatches where rnaSPAdes called something different than
// what is shown by a pileup of the reads themselves

process refineContigs {

    // Save consensus.fasta and binned-reads.fastq with their sample name & TVV species
    publishDir path: "${params.outputDirectory}/analysis/05_refinement/",
               pattern: "consensus.fasta.gz",
               mode: "copy",
               saveAs: { filename -> "${reads.getSimpleName()}.${filename}" }

   publishDir path: "${params.outputDirectory}/analysis/05_refinement/",
              pattern: "*fq",
              mode: "copy"

    input:
    tuple file(mapped_bam), file(contigs), file(reads) from mapped_bam

    output:
    tuple file("consensus.fasta.gz"), file(reads) into refined_contigs
    file reads

    """
    # Convert the TVV-aligned-reads (BAM) into a pileup (VCF)
    bcftools mpileup \
        -d 1000000 \
        -f $contigs \
        $mapped_bam > \
        pileup.vcf

    # Call the variants (BCF file)
    bcftools call -m -Ob -o variants_called.bcf pileup.vcf

    # Index the calls.bcf file
    bcftools index variants_called.bcf

    # Combine the reference fasta and the called-variants into a consensus FASTA
    bcftools consensus \
        -f $contigs \
        variants_called.bcf > \
        consensus.fasta

    # Compress the consensus fasta
    gzip consensus.fasta
    """

}

process classification {

    // Diamond parameters
    // block_size_to_use = ${available_memory}/10

    // Save classifications files
    publishDir path: "${params.outputDirectory}/analysis/06_classification/",
               pattern: "classification.txt",
               mode: "copy",
               saveAs: { filename -> "${reads.getSimpleName()}.${filename}" }

    // Take in refined contigs and reads files
    input:
    tuple file(refined_contigs), file(reads) from refined_contigs

    output:
    tuple file("classification.txt"), file(reads) into classified_contigs

    """
    # Run diamond
    diamond \
    blastx \
    --verbose \
    --more-sensitive \
    --db $params.diamondDB \
    --query $refined_contigs \
    --out classification.txt \
    --outfmt 102 \
    --max-hsps 1 \
    --top 1 \
    --block-size $params.blockSize \
    --index-chunks 2 \
    --threads $params.threads
    """

}

process taxonomy {

    // Save translated classification files containing the full taxonomic lineages
    publishDir path: "${params.outputDirectory}/analysis/06_classification/",
               pattern: "classification.taxonomy.txt",
               mode: "copy",
               saveAs: { filename -> "${reads.getSimpleName()}.${filename}" }

    input:
    tuple file(classifications), file(reads) from classified_contigs

    output:
    file "classification.taxonomy.txt"

    """
    diamondToTaxonomy.py $classifications
    """

}
