#!/usr/bin/env nextflow

// Setup default parameters if not provided any
//params.sample = "sample"
params.reads = "data/*_R{1,2}.fastq"
params.outputDirectory = "output/"
params.phiX = "genomes/phiX.fasta"
params.tvvDirectory = "genomes/"
params.diamondDB = "/n/data1/hms/mbib/nibert/austin/diamond/ncbi-viruses-bacteria-tv.dmnd"
params.blockSize = "2"
params.threads = "4"

// Read in the sample's reads files
paired_reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true, flat: true)

// Trim adapters from the reads
process trimAdapters {
    publishDir "${params.outputDirectory}/analysis/01_adapter_trimmed/", mode: "copy"

    input:
    tuple val(sampleID), file(forward_reads), file(reverse_reads) from paired_reads_ch

    output:
    tuple val(sampleID), file("${sampleID}*val_1.fq.gz"), file("${sampleID}*val_2.fq.gz") into trimmed_reads

    """
    bash trimAdapters.sh \
    -s $sampleID -1 $forward_reads -2 $reverse_reads -o "./"
    """
}

// Remove any reads that map to phiX (spiked in during RNA-seq library prep for sequencing
process depletePhiX {
    publishDir "${params.outputDirectory}/analysis/02_phiX_depleted/", mode: "copy"

    input:
    tuple val(sampleID), file(trimmed_fwd_reads), file(trimmed_rev_reads) from trimmed_reads

    output:
    tuple val(sampleID), \
          file("${sampleID}*R1.fq.gz"), \
          file("${sampleID}*R2.fq.gz") \
    into phiX_depleted_reads

    """
    bash depletePhiX.sh \
    -s $sampleID -1 $trimmed_fwd_reads -2 $trimmed_rev_reads -o "./" -p "${workflow.launchDir}/$params.phiX"
    """

}

process mapToTVV {
    publishDir "${params.outputDirectory}/analysis/03_binned_reads/fastq/", mode: "copy"

    input:
    tuple val(sampleID), file(forward_reads), file(reverse_reads) from phiX_depleted_reads

    output:
    file "${sampleID}_*fq" into binned_tvv, binned_tvv_fastq

    """
    bash mapToTVV.sh \
    -s $sampleID -1 $forward_reads -2 $reverse_reads  -o "./" -t "$workflow.launchDir/$params.tvvDirectory"
    """
}

// Wrangle all six (one per TVV species + satellites) files from binned_tvv & binned_tvv_fastq
binned_tvv = binned_tvv.flatten()
             .map { file -> tuple(file.simpleName, file)}

binned_tvv_fastq = binned_tvv_fastq.flatten()
                   .map { file -> tuple(file.simpleName, file)}


process fastqToFasta {
    publishDir "${params.outputDirectory}/analysis/03_binned_reads/fasta/", mode: "copy"

    input:
    tuple val(sample_and_tvv_species), file(fastq) from binned_tvv_fastq

    output:
    file "${sample_and_tvv_species}.fasta"

    """
    seqtk seq -A $fastq > "${sample_and_tvv_species}.fasta"
    """
}

process deNovoAssembly {

    // Now that we're at the multiple TVV files per sample stage, I am going to run these scripts
    // within their own sample+species specific work directories to prevent overwriting/interference,
    // but save them all together in the same directory at the end

    publishDir path: "${params.outputDirectory}/analysis/04_denovoassembly",
               pattern: "${sample_and_tvv_species}.transcripts.fasta",
               mode: "copy"

    input:
    tuple val(sample_and_tvv_species), file(reads) from binned_tvv

    output:
    tuple val(sample_and_tvv_species), \
          file("${sample_and_tvv_species}.transcripts.fasta"), \
          file(reads) \
    into tvv_contigs

    """
    # Build contigs with rnaSPAdes & drop any short contigs <300 nt;
    # if no contigs can be built that TVV species, output an empty file

    rnaspades.py -s $reads -o "${sample_and_tvv_species}/" && \
    seqtk seq -L 300 "${sample_and_tvv_species}/transcripts.fasta" > \
                     "${sample_and_tvv_species}/transcripts.trimmed.fasta" && \
    mv "${sample_and_tvv_species}/transcripts.trimmed.fasta" \
       "${sample_and_tvv_species}/transcripts.fasta" || \
    touch "${sample_and_tvv_species}/transcripts.fasta"

    # If rnaSPAdes partially succeeds (can only make low-quality contigs), it will end silently without providing a transcripts file
    if [[ ! -f "./${sample_and_tvv_species}/transcripts.fasta" ]]; then
        mkdir -p "./${sample_and_tvv_species}/"
        touch "./${sample_and_tvv_species}/transcripts.fasta"
    fi

    # Rename transcripts file to include sample name & move into main folder
    mv "${sample_and_tvv_species}/transcripts.fasta" "${sample_and_tvv_species}.transcripts.fasta"
    """
}

process refineContigs {

    // Map the TVV-reads to the denovo-assembled-contigs with BWA;
    // with the reads mapped to their denovo-assemblies, refine the assembly by
    // finding any mismatches where rnaSPAdes called something different than
    // what is shown by a pileup of the reads themselves

    // Save the refined TVV contigs
    publishDir path: "${params.outputDirectory}/analysis/05_refinement/",
               pattern: "${sample_and_tvv_species}.refined_contigs.fasta.gz",
               mode: "copy"

   // Save the variants-called bcf file
   publishDir path: "${params.outputDirectory}/analysis/05_refinement/mapping",
              pattern: "${sample_and_tvv_species}.variants_called.bcf",
              mode: "copy"

    // Save the BAM file with the name of the TVV species
    publishDir path: "${params.outputDirectory}/analysis/05_refinement/mapping/",
               pattern: "${sample_and_tvv_species}.reads_mapped_to_contigs.sorted.bam",
               mode: "copy"

    // Save the mapping-statistics file with the name of the TVV species
    publishDir path: "${params.outputDirectory}/analysis/05_refinement/mapping/",
              pattern: "${sample_and_tvv_species}.reads_mapped_to_contigs.sorted.stats",
              mode: "copy"

    // Only read in the files for TVV species where rnaSPAdes could actually construct contigs
    input:
    tuple val(sample_and_tvv_species), \
          file(contigs), \
          file(reads) \
    from tvv_contigs.filter { it.get(1).size() > 0 }

    output:
    tuple val(sample_and_tvv_species), \
          file("${sample_and_tvv_species}.refined_contigs.fasta.gz"), \
          file(reads) \
    into refined_contigs, refined_contigs_and_reads_for_coverage

    """
    bash refineContigs.sh  \
    -s "${sample_and_tvv_species}" -r $reads -c $contigs -o "./" -t $params.threads
    """
}

// Map the reads to the contigs to determine per-contig coverage
process coverage {

    publishDir path: "${params.outputDirectory}/analysis/06_coverage/",
               pattern: "${sample_and_tvv_species}.contigs_coverage.txt",
               mode: "copy"

    input:
    tuple val(sample_and_tvv_species), \
          file(refined_contigs), \
          file(reads) \
    from refined_contigs_and_reads_for_coverage

    output:
    tuple val(sample_and_tvv_species), \
          file(refined_contigs), \
          file("${sample_and_tvv_species}.contigs_coverage.txt") \
    into contigs_with_coverage

    """
    # Index contigs for BWA
    bwa index -p "${sample_and_tvv_species}_index" $refined_contigs

    # Map reads to contigs with BWA-mem
    bwa mem -t $params.threads "${sample_and_tvv_species}_index" $reads | \
    samtools sort --threads $params.threads -o "${sample_and_tvv_species}.mapped.bam"

    # Calculate the mean-depth (i.e., coverage) per contig; keep each contig's name & coverage; throw away header; sort by coverage
    samtools coverage "${sample_and_tvv_species}.mapped.bam" | \
    cut -f 1,7 > "${sample_and_tvv_species}.contigs_coverage.txt"
    """
}

process classification {

    // Save classifications files
    publishDir path: "${params.outputDirectory}/analysis/07_classification/",
               pattern: "${sample_and_tvv_species}.classification.txt",
               mode: "copy"

    // Take in refined contigs and reads files only if it's from the unmapped read/contigs
    input:
    tuple val(sample_and_tvv_species), \
          file(refined_contigs),  \
          file(coverage) \
    from contigs_with_coverage.filter { it.get(0) =~/unmapped/ }

    output:
    tuple val(sample_and_tvv_species), \
          file("${sample_and_tvv_species}.classification.txt"), \
          file(coverage) \
    into classified_contigs

    """
    # Run diamond
    diamond \
    blastx \
    --verbose \
    --more-sensitive \
    --db $params.diamondDB \
    --query $refined_contigs \
    --out "${sample_and_tvv_species}.classification.txt" \
    --outfmt 6 qseqid staxids evalue bitscore pident qcovhsp \
    --max-hsps 1 \
    --top 0 \
    --block-size $params.blockSize \
    --index-chunks 2 \
    --threads $params.threads

    """
}

process taxonomy {

    // Save translated classification files containing the full taxonomic lineages
    publishDir path: "${params.outputDirectory}/analysis/08_taxonomy/",
               pattern: "${sample_and_tvv_species}.classification.taxonomy.txt",
               mode: "copy"

   publishDir path: "${params.outputDirectory}/analysis/08_taxonomy/",
              pattern: "${sample_and_tvv_species}.final_table.txt",
              mode: "copy"

    input:
    tuple val(sample_and_tvv_species), file(classifications), file(coverage) from classified_contigs

    output:
    file "${sample_and_tvv_species}.final_table.txt"

    """
    # Translate the DIAMOND results to full lineages
    diamondToTaxonomy.py $classifications

    # Join the coverage values and the taxonomy results
    join \
        -j 1 \
        -t \$'\t' \
        --check-order \
        <(sort -k1,1 $coverage) \
        <(grep -v "^#" "${sample_and_tvv_species}.classification.taxonomy.txt" | sort -k1,1) | \
    sort -rgk2,2 > \
    "${sample_and_tvv_species}.contigs_coverage_taxonomy.txt"

    # Make a header for a final results table
    echo -e \
        "#Contig\t" \
        "#Coverage\t" \
        "#TaxonID\t" \
        "#e-value\t" \
        "#Bitscore\t" \
        "#PercentIdentity\t" \
        "#QueryCoverage\t" \
        "#Domain\t" \
        "#Kingdom\t" \
        "#Phylum\t" \
        "#Class\t" \
        "#Order\t" \
        "#Family\t" \
        "#Genus_species" \
    > "${sample_and_tvv_species}.final_table.txt"

    # Add the data to the final with just the header
    cat "${sample_and_tvv_species}.contigs_coverage_taxonomy.txt" >> \
        "${sample_and_tvv_species}.final_table.txt"
    """

}
