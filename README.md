# Fasta-Region-Inspector

## Introduction
Fasta-Region-Inspector (**FRI**) is a bioinformatics tool for analyzing somatic hypermutation (SHM).

## Purpose
Analyzing cancer cohorts for SHM is a interesting cancer-genomics related research question, and one that is computationally intensive.  SHM occurs when immunoglobulin genes accumulate apparently random point mutations within productively rearranged V, D, and J segments, as defined by [Elaine S. Jaffe MD, in Hematopathology, 2017](https://www.sciencedirect.com/topics/immunology-and-microbiology/somatic-hypermutation).  More specifically, [Blood. 2009 Apr 16; 113(16): 3706–3715](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2670789/) states that the hallmark of SHM is the increased percentage of mutations in hypermutable RGYW (underlined residue is the most frequently mutated) motifs.  Correctly identifying SHM in large cancer genomics datasets can be a complex and tedious computational problem. 

This bioinformatics tool gives researchers the ability to answer the following questions:
- Which variants are within 2 kb of the transcription start site (TSS) of the corresponding gene?
- Where do user-defined mapped ambiguity strings lie within 2 kb of the TSS?
- Amalgamating the previous 2 questions, which variant(s) are found within a mapped ambiguity string, that also lie within 2 kb of the TSS?

This tool aims to answer common SHM variant-level questions in software package that provides:
- Excellent runtime performance in a robust, functional implementation.
- Minimized memory usage.
- A simple, YAML input file format.
- Clean, informative stdout logging.

## Theory and Implementation

## Improvements from Previous Implemenation
This version of **FRI** is vastly improved upon from the [old version](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD) in numerous ways, including:
- A completely re-worked internal memory representation of user-defined fasta data.
  - The old implementation loaded the fasta file in a strict ByteString, and then parsed it into an appropriate data structure using the [elynx-seq](https://hackage.haskell.org/package/elynx-seq) library.
  - There's nothing inherently wrong with the approach, but it doesn't scale well with sequencing data as this type of data is typically very large.  That means the GC has to walk through this in-memory structure many times for all the functions that operate on this data.
    - These functions include:
      - [ambiguityCodesWithinRegionCheck](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L381)
      - [ambiguityCodesWithinRegionCheckSmall](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L396)
      - [subStrLocations](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L458)
      - [subStrLocationsSmallForward](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L487)
      - [subStrLocationsSmallReverse](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L471)
      - [grabFastaSequence](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L508)
      - [smallGrabFastaSequence](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L514)
  - Previously, the [old version](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD) consumed memory in such a way that (using a concatenated amalgamation of all GRCh38 homo_sapiens genome assembly chromosome-level fastas) during a typical run that you would be required to have access to and run the software with at least 30+ GB of memory.
- A new, YAML input file format.
  - The new [YAML](https://yaml.org/) input file format

