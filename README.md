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

## Improvements from Previous Implemenation
This version of **FRI** is vastly improved upon from the [old version](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD) in numerous ways, including:
- Migration from single script to [Stack](https://docs.haskellstack.org/en/stable/README/) project
  - Stack is a cross-platform program for developing Haskell projects.
  - It creates completely reproducible and robust development environments.
  - This makes ongoing maintenance, updates, bugfixes, any and all changes to the codebase a breeze.
- A completely re-worked internal memory representation of user-defined FASTA data using [compact regions](http://ezyang.com/papers/ezyang15-cnf.pdf).
  - The old implementation loaded the FASTA file in a strict ByteString, and then parsed it into an appropriate data structure using the [elynx-seq](https://hackage.haskell.org/package/elynx-seq) library.
  - There's nothing inherently wrong with the approach, but it doesn't scale well with sequencing data as this type of data is typically very large.  That means the garbarge collector (GC) has to walk through this in-memory structure many times for all the functions that operate on this data.
    - These functions include:
      - [ambiguityCodesWithinRegionCheck](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L381)
      - [ambiguityCodesWithinRegionCheckSmall](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L396)
      - [subStrLocations](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L458)
      - [subStrLocationsSmallForward](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L487)
      - [subStrLocationsSmallReverse](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L471)
      - [grabFastaSequence](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L508)
      - [smallGrabFastaSequence](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L514)
  - Previously, the [old version](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD) consumed memory in such a way that (using a concatenated amalgamation of all GRCh38 homo_sapiens genome assembly chromosome-level fastas) during a typical run that you would be required to have access to and run the software with at least 40+ GB of memory.
  - This same run now only requires about as much memory as the size of the FASTA file plus 1/2.
- A new, YAML input file format.
  - The new [YAML](https://yaml.org/) input file format replaces the [custom command-line argument string](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L840) required to run the [old version](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD) with a new simple key-value format:
    - ```
      Variants:
        - Sample: 'SMP1'
          HGNC_Symbol: 'AFF3'
          Chromosome: 'chr2'
          Start_Position: '100007433'
          End_Position: '100007433'
          Reference_Allele: 'G'
          Alternate_Allele: 'A'
          ENST: 'ENST00000409579'
        - Sample: 'SMP2'
        ...
      ```
  - The rest of the command-line arguments are replaced with appropriate compact nested mappings (see [Example 2.12 Compact Nested Mapping](https://yaml.org/spec/1.2.2/#chapter-2-language-overview))
- Switching from the [Copying Collector GC](https://gitlab.haskell.org/ghc/ghc/-/wikis/commentary/rts/storage/gc/copying) to the new [Non-Moving GC](https://www.cs.unh.edu/~dietz/papers/gamari2020alligatordemo.pdf)
  - [GHC](https://www.haskell.org/ghc/), the the de-facto compiler for the Haskell programming language.
  - GHC initially had only one implementation of GC, namely the Copy Collector GC.
    - This form of GC has a stop-the-world approach, cause the program to pause as the memory usage raises to a level where de-allocation is necessary.
  - GHC now has a new implemenation of the GC, the Non-Moving GC, which this program utilizes.
    - The [new implemenation](https://gitlab.haskell.org/ghc/ghc/-/commit/7f72b540288bbdb32a6750dd64b9d366501ed10c) of GC has a concurrent mark & sweep garbage collector to manage the old generation. The concurrent nature of this collector typically results in significantly reduced maximum and mean pause times in applications with large working sets.
      - For more details, please see the following [presentation](https://bgamari.github.io/media/2018-11-18-nonmoving-gc-for-ghc.pdf).
  - You can switch between the two forms of GC using the GHC command-line argument ```-with-rtsopts=-xn```.

## Configuration YAML

**FRI** utilizes a configuration YAML file to provide the necessary components for a successful run.

The following keys are **required**:
- ```Fasta``` -> The filepath to the FASTA file. (String)
- ```Variants``` -> A [compact nested mappings](https://yaml.org/spec/1.2.2/#chapter-2-language-overview) (see [below](https://github.com/Matthew-Mosior/Fasta-Region-Inspector/blob/main/README.md#variant-type) for more information).
- ```Ambiguity_Codes``` -> An array of ambiguity codes to search against. (String)
- ```Output_Directory``` -> The filepath to the output directory (must already exist). (String)
- ```Keep_BioMart``` -> Whether or not to keep the BioMart file from the wget system call. (Boolean)
- ```Ignore_Strandedness``` -> Whether or not to ignore strandness of the respective gene and search the TSS in both directions. (Boolean)

The following keys are **optional**:
- ```TSS_Window_Size``` -> The TSS window size to search across.

## Variant Type
