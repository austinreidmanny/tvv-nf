#!/usr/bin/env nextflow

/*
   ------------------------------------------------------------------------------------------------
   WELCOME TO THE tvv-nf PIPELINE
   ------------------------------------------------------------------------------------------------
   Objective:

   This pipeline will take in a set of paired-end dsRNA-seq reads, map them to all TVV species,
   produce refined contigs (attempting to generate full-length viral genomes), and will assemble all
   non-TVV reads into contigs to look for new viruses in these Trichomonas vaginalis samples
   ------------------------------------------------------------------------------------------------
   Software:

   Download the latest version of tvv-nf @ https://github.com/austinreidmanny/tvv-nf
   Dependencies & requirements are
     * Operating system: Linux/macOS
     * Nextflow (https://www.nextflow.io)
     * Conda (https://docs.conda.io/en/latest/miniconda.html)
     * DIAMOND database
     * Kraken database
   ------------------------------------------------------------------------------------------------
   Usage:

   nextflow run tvv_dsRNAseq_analysis.nf --reads "data/*R{1,2}.fastq.gz"

   acceptable '--reads' formatting:
     "*R{1,2}.fastq"
     "sample*R{1,2}.fastq"
     "directory/*R{1,2}.fastq"
     "directory/sample*R{1,2}.fastq"
     (any of the above plus ".gz")

   optional parameters:
   --outputDirectory "output/dir/"
   --diamondDB "path/to/db"
   --krakenDB "path/to/db"
   --blockSize "#" [for DIAMOND, should be approx 1/10 of your memory/RAM]
   --threads "#"
   --phiX "custom/path/to/phiX.fasta"
   --tvvDirectory "custom/path/to/tvv*.fasta"
   ------------------------------------------------------------------------------------------------
   Contact:

   Austin R. Manny
   Nibert Lab @ Harvard Medical School
   austinmanny@g.harvard.edu
   github.com/austinreidmanny
   ------------------------------------------------------------------------------------------------
*/

//================================================================================================//

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
    into phiX_depleted_reads, cleaned_reads_for_classification

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
    file "${sampleID}_*fq" into binned_tvv_for_fasta_conversion
    file "${sampleID}_*R1.fq" into binned_forward_reads
    file "${sampleID}_*R2.fq" into binned_reverse_reads

    """
    bash mapToTVV.sh \
    -s $sampleID -1 $forward_reads -2 $reverse_reads  -o "./" -t "$workflow.launchDir/$params.tvvDirectory"
    """
}


// ---------------------------------------------------------------------------------------------- //
// Wrangle all six (one per TVV species + satellites)
// ---------------------------------------------------------------------------------------------- //

// Identify the TVV viral each file is binned to
binned_forward_reads = binned_forward_reads.flatten()
                       .map { file -> tuple(file.simpleName, file) }
binned_reverse_reads = binned_reverse_reads.flatten()
                       .map { file -> tuple(file.simpleName, file) }

// Match the corresponding forward and reverse reads
(binned_tvv, binned_tvv_for_assembly) = binned_forward_reads.combine(binned_reverse_reads, by: 0).into(2)

// For FASTA conversion, no pair matching is required
binned_tvv_for_fasta_conversion = binned_tvv_for_fasta_conversion.flatten()
                                  .map { file -> tuple(file.simpleName, file)}

// ---------------------------------------------------------------------------------------------- //

process fastqToFasta {
    publishDir "${params.outputDirectory}/analysis/03_binned_reads/fasta/", mode: "copy"

    input:
    tuple val(sample_and_tvv_species), file(fastq) from binned_tvv_for_fasta_conversion

    output:
    file "*.fasta"

    """
    seqtk seq -A $fastq > "${fastq.baseName}.fasta"
    """
}

// Indicate forward or reverse read in read headers by changing spaces to underscores
process fixHeaders {
    publishDir "${params.outputDirectory}/analysis/04_cleanup_binned_reads/01_fix_headers/", mode: "copy"

    input:
	tuple val(sample_and_tvv_species), file(forward_reads), file(reverse_reads) from binned_tvv

    output:
    tuple val(sample_and_tvv_species), file("${sample_and_tvv_species}.fixed-headers_R1.fasta.gz"), file("${sample_and_tvv_species}.fixed-headers_R2.fasta.gz") into fixed_reads

    shell:
    $/
    seqtk seq -A $forward_reads | sed 's/^>/>R1_/g' | gzip -c > "${sample_and_tvv_species}.fixed-headers_R1.fasta.gz"
    seqtk seq -A $reverse_reads | sed 's/^>/>R2_/g' | gzip -c > "${sample_and_tvv_species}.fixed-headers_R2.fasta.gz"
    /$

}

process combineFasta {
	publishDir "${params.outputDirectory}/analysis/04_cleanup_binned_reads/02_combine_reads/", mode: "copy"

	input:
	tuple val(sample_and_tvv_species), file(forward_reads), file(reverse_reads) from fixed_reads

	output:
	tuple val(sample_and_tvv_species), file("*.combined.fasta") into combined_fasta

	"""
	zcat $forward_reads $reverse_reads > "${sample_and_tvv_species}.combined.fasta"
	"""
}

process cleanupBinnedReads {
	publishDir "${params.outputDirectory}/analysis/04_cleanup_binned_reads/03_cleaned_up/", mode: "copy"

    input:
  	tuple val(sample_and_tvv_species), file(reads) from combined_fasta.filter{ name, file -> name =~ /tvv\d/ } // ignore unmapped reads

  	output:
  	file "*summary.txt"
  	file "*results.txt" into cleaned_results

  	"""
  	bash cleanup_blast.sh \
  	-i $reads -p "${workflow.launchDir}/resources/cleanup_blast/minimal_tvv_db" -o "./" -m 8 -e 1e-3
	   """
}



process parseCleanupResults {
    publishDir "${params.outputDirectory}/analysis/04_cleanup_binned_reads/04_parsed/", mode: "copy"

    input:
    file results from cleaned_results

    output:
    file "*parsed.txt"

    """
    bash parse_cleanup_results.sh \
	-f $results
    """

}

process deNovoAssembly {

    // Now that we're at the multiple TVV files per sample stage, I am going to run these scripts
    // within their own sample+species specific work directories to prevent overwriting/interference,
    // but save them all together in the same directory at the end

    publishDir path: "${params.outputDirectory}/analysis/05_denovoassembly",
               pattern: "${sample_and_tvv_species}.transcripts.fasta",
               mode: "copy"

    input:
    tuple val(sample_and_tvv_species), file(forward_reads), file(reverse_reads) from binned_tvv_for_assembly
        .filter{ name, fwd_reads_file, rev_reads_file -> name =~ /unmapped/ } // only assemble non-TVV reads

    output:
    tuple val(sample_and_tvv_species), file("${sample_and_tvv_species}.transcripts.fasta") into tvv_contigs

    """
    # --------------------------------------------------------------------------
    # Build contigs with rnaSPAdes & drop any short contigs <300 nt;
    # if no contigs can be built that TVV species, output an empty file
    # --------------------------------------------------------------------------

    # SPAdes does not play well with caching, will terminate if any files are in the output directory, so if it exists, just delete it
    rm -R "spades-output/" || true

    # Run SPAdes
    rnaspades.py -1 $forward_reads -2 $reverse_reads -o "spades-output/" --threads $task.cpus && \
        seqtk seq -L 300 "spades-output/transcripts.fasta" > "transcripts.trimmed.fasta" && \
        mv "transcripts.trimmed.fasta" "transcripts.fasta" || \
    touch "transcripts.fasta"

    # Rename transcripts file to include sample name
    mv "transcripts.fasta" "${sample_and_tvv_species}.transcripts.fasta"
    """

}

process classifyContigs {

    // Save classifications files
    publishDir path: "${params.outputDirectory}/analysis/06_classify_assemblies/",
               pattern: "${sample_and_tvv_species}.classification*txt",
               mode: "copy"

    // Take in refined contigs and reads files only if it's from the unmapped read/contigs
    input:
    tuple val(sample_and_tvv_species), file(contigs) from tvv_contigs

    output:
    file "${sample_and_tvv_species}.classification.txt"
    file "${sample_and_tvv_species}.classification.taxonomy.txt"

    """
    # --------------------------------------------------------------------------
    # Run diamond
    # --------------------------------------------------------------------------
    diamond \
    blastx \
    --verbose \
    --more-sensitive \
    --db $params.diamondDB \
    --query $contigs \
    --out "${sample_and_tvv_species}.classification.txt" \
    --outfmt 102 \
    --top 1 \
    --block-size $params.blockSize \
    --index-chunks 2 \
    --threads $task.cpus

    # --------------------------------------------------------------------------
    # Convert taxonomy IDs to useful lineages
    # --------------------------------------------------------------------------
    diamondToTaxonomy.py "${sample_and_tvv_species}.classification.txt"
    """
}

process classifyReads {
    // Now that the contigs are assembled and classified, I would like to also do a metatranscriptomic
    // census of just the (lightly processed/cleaned) unassembled reads

    publishDir path: "${params.outputDirectory}/analysis/07_classify_reads/",
               pattern: "${sampleID}.kraken-report.txt",
               mode: "copy"

    input:
    tuple val(sampleID), file(forward_reads), file(reverse_reads) from cleaned_reads_for_classification

    output:
    file "${sampleID}.kraken-report.txt"

    """
    kraken2 \
    --db $params.krakenDB \
    --paired --gzip-compressed --memory-mapping \
    --threads $task.cpus \
    --output ${sampleID}.kraken-output.txt \
    --report ${sampleID}.kraken-report.txt \
    $forward_reads \
    $reverse_reads
    """

}
