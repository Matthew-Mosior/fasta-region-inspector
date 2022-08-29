# Fasta-Region-Inspector

## Introduction
Fasta-Region-Inspector (**FRI**) is a bioinformatics tool for analyzing somatic hypermutation (SHM).

## Purpose
Analyzing cancer cohorts for SHM is a interesting cancer-genomics related research question, and one that is computationally intensive.  SHM occurs when immunoglobulin genes accumulate apparently random point mutations within productively rearranged V, D, and J segments, as defined by [Elaine S. Jaffe MD, in Hematopathology, 2017](https://www.sciencedirect.com/topics/immunology-and-microbiology/somatic-hypermutation).  More specifically, [Blood. 2009 Apr 16; 113(16): 3706–3715](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2670789/) states that the hallmark of SHM is the increased percentage of mutations in hypermutable RGYW motifs.  Correctly identifying SHM in large cancer genomics datasets can be a complex and tedious computational problem. 

This bioinformatics tool gives researchers the ability to answer the following questions:
- Which variants are within 2 kb of the transcription start site (TSS) of the corresponding gene?
- Where do user-defined mapped ambiguity strings lie within 2 kb of the TSS?
- Amalgamating the previous 2 questions, which variant(s) are found within a mapped ambiguity string, that also lie within 2 kb of the TSS?

## Improvements from Previous Implemenation

### Overview
This tool aims to answer common SHM variant-level questions in a software package that provides:
- Excellent runtime performance in a robust, functional implementation.
- Minimized memory usage.
- A simple, YAML input file format.
- Clean, informative stdout logging.

### Old vs. New
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
  - This same run now only requires about as much memory as ~2X the size of the FASTA file.
  - The minimized memory requirement is based on the usage of [compact regions](http://ezyang.com/papers/ezyang15-cnf.pdf).
    - [compact regions](http://ezyang.com/papers/ezyang15-cnf.pdf) serve two purposes (per [Data.Compact](https://hackage.haskell.org/package/compact-0.2.0.0/docs/Data-Compact.html) documentation):
      - Data stored in a Compact has no garbage collection overhead. The garbage collector considers the whole Compact to be alive if there is a reference to any object within it.
      - A Compact can be serialized, stored, and deserialized again. The serialized data can only be deserialized by the exact binary that created it, but it can be stored indefinitely before deserialization.
- A new, YAML input file format.
  - The new [YAML](https://yaml.org/) input file format replaces the [custom command-line argument string](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD/blob/master/src/fri.hs#L840) required to run the [old version](https://github.com/Matthew-Mosior/Fasta-Region-Inspector-OLD) with a simple key-value format:
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
  - This also supports a more legible, reproducible methodology since you have an easy to read and run file format.
- Switching from the [Copying Collector GC](https://gitlab.haskell.org/ghc/ghc/-/wikis/commentary/rts/storage/gc/copying) to the new [Non-Moving GC](https://www.cs.unh.edu/~dietz/papers/gamari2020alligatordemo.pdf)
  - [GHC](https://www.haskell.org/ghc/) is the the de-facto compiler for the Haskell programming language.
  - GHC initially had only one implementation of GC, namely the Copy Collector GC.
    - This form of GC has a stop-the-world approach, which causes the program to pause as the memory usage raises to a level where de-allocation is necessary.
  - GHC now has a new implemenation of the GC, the Non-Moving GC, which this program utilizes **by default**.
    - The [new implemenation](https://gitlab.haskell.org/ghc/ghc/-/commit/7f72b540288bbdb32a6750dd64b9d366501ed10c) of GC has a concurrent mark & sweep garbage collector to manage the old generation. The concurrent nature of this collector typically results in significantly reduced maximum and mean pause times in applications with large working sets.
      - For more details, please see the following [presentation](https://bgamari.github.io/media/2018-11-18-nonmoving-gc-for-ghc.pdf).
  - You can switch between the two forms of GC using the GHC command-line argument ```-with-rtsopts=-xn```.
- Dramatic runtime performance improvement
  - After all of the above changes (moving to stack or the new YAML input file doesn't affect runtime performance, but the rest of the points do), there is ~160X speedup in runtime (*benchmarks up-coming*) in the examples run locally so far (comparable tests are showing ~8 hrs down to ~3 minutes).
  - In general, the more ambiguity codes and variants there are to be analyzed by **FRI**, the larger the runtime disparity between the old and new version.

### Conclusion
Don't use the old version!

## Configuration YAML

**FRI** utilizes a configuration YAML file to provide the necessary components for a successful run.

The following keys are **required**:
- ```Fasta``` -> The filepath to the FASTA file. (String)
  - Please prefix each chromosome in the FASTA file with **chr**.
    - Ex. Chromosome 1 should look like -> chr1
- ```Variants``` -> A [compact nested mapping](https://yaml.org/spec/1.2.2/#chapter-2-language-overview) (see [below](https://github.com/Matthew-Mosior/Fasta-Region-Inspector/blob/main/README.md#variant-type) for more information).
- ```Ambiguity_Codes``` -> An array of ambiguity codes to search against. (String)
- ```Output_Directory``` -> The filepath to the output directory (must already exist). (String)
- ```Keep_BioMart``` -> Whether or not to keep the BioMart file from the wget system call. (Boolean)
- ```Ignore_Strandedness``` -> Whether or not to ignore strandness of the respective gene and search the TSS in both directions. (Boolean)
- ```Write_Ambiguity_Codes``` -> Whether or not to write CSV file containing all locations of user-defined mapped ambiguity code strings in TSS window. (Boolean)

The following keys are **optional**:
- ```TSS_Window_Size``` -> The TSS window size to search across. (String)

## Variant Type
The ```Variants``` [compact nested mapping](https://yaml.org/spec/1.2.2/#chapter-2-language-overview) key represents the required information surrounding a variant-of-interest to be examined by **FRI**.

The following keys are **required**:
- ```Sample``` -> The associated sample identifier for the respective variant (String).
- ```HGNC_Symbol``` -> The HGNC symbol for the respective variant (String).
- ```Chromosome``` -> The chromosome the respective variant lies within (String).
- ```Start_Position``` -> The start position for the respective variant (String).
- ```End_Position``` -> The end position for the respective variant (String).
- ```Reference_Allele``` -> The reference allele for the respective variant (String).
- ```Alternate_Allele``` -> The alternate allele for the respective variant (String).
- ```ENST``` -> The ENST for the respective variant (String).

### Notes
The above data can typically be easily grepped/programmed for from the output of a bioinformatics pipeline run, such as the [alignment_exome.cwl](https://github.com/genome/analysis-workflows/blob/master/definitions/pipelines/alignment_exome.cwl).

You can also start from a [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)-annotated variant file, which can typically be created once you have run an aligner/variant caller (see above).  Starting from a VEP-annotated file will allow you to filter the variants on gene so that you can better target your analysis with this tool.

## Building the project
This software was developed on an M1 mac.

It has not been tested on other operating systems or chipsets, but should run perfectly fine on [Stack](https://docs.haskellstack.org/en/stable/README/)-supported operating systems, including:
- Windows
- MacOS
- Ubuntu
- Debian
- Fedora
- Arch Linux
- CentOS
- Red Hat
- Amazon Linux
- openSUSE
- NixOS
- Generic Linux


If you are on an M1 mac, please build the project using the following:
```
% stack build --arch aarch64
```

If you are **NOT** on an M1 mac, please build the project using the following:
```
% stack build
```

## Example Usage
**FRI** is easy to use, as it only requires a single command-line positional argument, the configuration YAML.

If you are on an M1 mac, please run the project using the following:
```
% stack exec --arch aarch64 fasta-region-inspector-exe /path/to/configuration.yaml
```

If you are **NOT** an M1 mac, please run the project using the following:
```
% stack exec fasta-region-inspector-exe /path/to/configuration.yaml
```

## Dependency - Theorem Solvers
The program uses the [sbv](https://hackage.haskell.org/package/sbv) package, which interfaces with theorem solver(s).

This package interfaces with one of following supported theorem solvers:
- [ABC](http://www.eecs.berkeley.edu/~alanmi/abc/) from University of Berkeley.
- [Boolector](http://fmv.jku.at/boolector/) from Johannes Kepler University.
- [Bitwuzla](http://bitwuzla.github.io/) from Stanford University.
- [CVC4](http://cvc4.github.io/) and [CVC5](http://cvc5.github.io/) from Stanford University and the University of Iowa.
- [DReal](http://dreal.github.com/) from CMU.
- [MathSAT](http://mathsat.fbk.eu/) from FBK and DISI-University of Trento.
- [Yices](http://yices.csl.sri.com/) from SRI.
- [Z3](http://github.com/Z3Prover/z3/wiki) from Microsoft.

This software has been tested with both Z3 and Yices, although the other theorem solvers listed above should work just as well.

You will need to install and set-up one of the above theorem solvers in order for **FRI** to calculate the mapped ambiguity code string for the given ambiguity code(s).

Please see http://leventerkok.github.io/sbv/ to get more information on how theorem solvers work.

## Output
**By default**, **FRI** outputs a single CSV file,namely *variants_in_ambiguity_codes.csv*, containing all of the information needed to answer the questions posed at the top of this README.

Each user-defined variant will be included in this CSV file, with the following fields:
- ```Variant``` -> A description of the user-defined variant.
- ```Region``` -> A description of the region related metadata pulled from BioMart for the user-defined variant.
- ```Variant_Within_Region``` -> A boolean (Y/N) of whether or not the user-defined variant lies within 2 kb of the TSS.
- ```Ambiguity_Code``` -> The user-defined ambiguity code used for the stringsearch.
  - This can be the reverse complement of the user-defined ambiguity code, depending on the strand directionality from the BioMart for the variant.
- ```Mapped_Nucleotide_String``` -> The [Satisfiability Modulo Theories](https://hackage.haskell.org/package/sbv) calculated mapped ambuguity code string for the given ambiguity code.
- ```Ambiguity_Code_String_Locations_Within_TSS``` -> The starting position of the first base pair for the mapped nucleotide string.
  - Again, the strand directionality can change which way the mapped nucleotide string is interpreted (see below).

### Notes

```N/A``` values in the ```Ambiguity_Code```, ```Mapped_Nucleotide_String```, or ```Ambiguity_Code_String_Locations_Within_TSS``` fields indicate that the respective variant did not lie within a mapped nucleotide string within 2 kb of the TSS.

Depending on how you set the ```Write_Ambiguity_Codes``` key in the configuration YAML, you can also get an additional CSV file containing all locations of user-defined mapped ambiguity code strings in the TSS window.

## Visualizing the Output
Interpretation of the *variants_in_ambiguity_codes.csv* file is a simple process.

The following example will help with ground your understanding of the output.

### Note

If you wish to follow along with your own data, please download [IGV](https://software.broadinstitute.org/software/igv/).

IGV is a helpful resource for visualizing your results (all of the following screenshots in this section are from IGV). 

### Examine the Mapped Nucleotide String<sup>**Assumes a working knowledge of IGV**</sup>

|Variant|Region|Variant_Within_Region|Ambiguity_Code|Mapped_Nucleotide_String|Ambiguity_Code_String_Locations_Within_TSS|
|-------|------|---------------------|--------------|------------------------|------------------------------------------|
|SMP1:H1-5:chr6:27867219:27867219:G:A:ENST00000331442|6:27867588:&#x200B;-1:H1-5|Y|GYW|GCT|27867220|

Let's take a look at what this particular variant looks like in IGV:
<br></br>
<p align="center">
  <img src="/img/forwardstrand.png">
</p>
<p align="center">
  The mapped nucleotide string on the forward strand.
</p>
<p align="center">
  <img src="/img/reversestrand.png" />
</p>
<p align="center">
  The mapped nucleotide string on the reverse strand.
</p>
 
In IGV, you are able to view a nucleotide sequence in both the forward strand orientation (**1**) and the reverse strand orientation (**-1**) (to do this, click on the arrow icon on the bottom left portion of IGV near **Sequence**).

You can see in the example above that the variant is found on the **reverse** strand.

Because of this, you should focus on the reverse strand within IGV (see the reverse strand image above).

To help orient yourself, the **T** is at position **27867218** and the **G** is at position **27867220** on chromosome 6. 

Our variant-of-interest is the **C**, at position **27867219**, due to it being on the reverse strand (the complement to **G** is **C**, and the input to the tool's ```Reference_Allele``` and ```Alternate_Allele``` are in reference to the forward strand).

Now that you have a better understanding of what the output means, you can be confident in filtering the output file in a way that is congruent with the research question you are trying to answer.

## Example Stdout (logging)
The following is a real stdout (log) of a **FRI** run:

```
        ______           __           ____             _                ____                           __
       / ____/___ ______/ /_____ _   / __ \___  ____ _(_)___  ____     /  _/___  _________  ___  _____/ /_____  _____
      / /_  / __ `/ ___/ __/ __ `/  / /_/ / _ \/ __ `/ / __ \/ __ \    / // __ \/ ___/ __ \/ _ \/ ___/ __/ __ \/ ___/
     / __/ / /_/ (__  ) /_/ /_/ /  / _, _/  __/ /_/ / / /_/ / / / /  _/ // / / (__  ) /_/ /  __/ /__/ /_/ /_/ / /
    /_/    \__,_/____/\__/\__,_/  /_/ |_|\___/\__, /_/\____/_/_/_/__/___/_/ /_/____/ .___/\___/\___/\__/\____/_/
                                             /____/_/ __ \ <  // __ \ / __ \      /_/
                                             | | / / / / / / // / / // / / /
                                             | |/ / /_/ / / // /_/ // /_/ /
                                             |___/\____(_)_(_)____(_)____/

                                           Copyright (c) Matthew C. Mosior 2022

[2022-08-18 14:05:43.530389 EDT] Starting up Fasta Region Inspector v0.1.0.0 ...
[2022-08-18 14:05:43.531635 EDT] Query BioMart for regions data ...
[2022-08-18 14:05:43.531702 EDT] Generating BioMart compatible XML ...
[2022-08-18 14:05:43.594794 EDT] Querying and downloading region data from BioMart via HTTP request ...
[2022-08-18 14:05:52.646374 EDT] Successfully queried and returned region data from BioMart via HTTP request ...
[2022-08-18 14:05:52.646695 EDT] Writing BioMart region data to file biomartresult.txt in output directory ...
[2022-08-18 14:05:52.647817 EDT] Determining whether each variant is within the its respective genes TSS ...
[2022-08-18 14:05:52.647948 EDT] Calculating the reverse complement of each user defined ambiguity code ...
[2022-08-18 14:05:52.648125 EDT] Creating list of tuples to define directionality of each forward strand ambiguity code ...
[2022-08-18 14:05:52.648235 EDT] Creating list of tuples to define directionality of each reverse strand ambiguity code ...
[2022-08-18 14:05:52.648295 EDT] Generating all possible ambiguity codes strings using SMT solver ...
[2022-08-18 14:05:52.734363 EDT] Preparing ambiguity code strings to determine whether each lies within its respective TSS ...
[2022-08-18 14:05:52.734412 EDT] Reading fasta file into strict ByteString ...
[2022-08-18 14:05:55.160107 EDT] Parsing fasta file into [Sequence] ...
[2022-08-18 14:06:48.863554 EDT] Extract every sequences name and characters from cfasta ...
[2022-08-18 14:06:48.863696 EDT] Putting namesandcharacters into compact region ...
[2022-08-18 14:08:47.505148 EDT] Determing whether each ambiguity code string lies within its respective TSS ...
[2022-08-18 14:08:47.530243 EDT] Processing region data associated with gene BTG1 ...
[2022-08-18 14:08:47.530264 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and BTG1 strand orientation is -1 ...
[2022-08-18 14:08:47.533043 EDT] Processing region data associated with gene MYH10 ...
[2022-08-18 14:08:47.533057 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and MYH10 strand orientation is -1 ...
[2022-08-18 14:08:47.533085 EDT] Processing region data associated with gene ARID5B ...
[2022-08-18 14:08:47.533098 EDT] Extracting fasta sequence associated with chromosome: 10, tss: 61901699, strand: 1, gene: ARID5B ...
[2022-08-18 14:08:47.533115 EDT] Processing mapped ambiguity code TGC ...
[2022-08-18 14:08:47.535758 EDT] Processing mapped ambiguity code TAC ...
[2022-08-18 14:08:47.535785 EDT] Processing mapped ambiguity code AGC ...
[2022-08-18 14:08:47.542388 EDT] Processing mapped ambiguity code AAC ...
[2022-08-18 14:08:47.542420 EDT] Processing region data associated with gene EBF1 ...
[2022-08-18 14:08:47.542447 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and EBF1 strand orientation is -1 ...
[2022-08-18 14:08:47.542482 EDT] Processing region data associated with gene H1-5 ...
[2022-08-18 14:08:47.542496 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and H1-5 strand orientation is -1 ...
[2022-08-18 14:08:47.542518 EDT] Processing region data associated with gene H1-2 ...
[2022-08-18 14:08:47.542528 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and H1-2 strand orientation is -1 ...
[2022-08-18 14:08:47.542549 EDT] Processing region data associated with gene SGK1 ...
[2022-08-18 14:08:47.542559 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and SGK1 strand orientation is -1 ...
[2022-08-18 14:08:47.542630 EDT] Processing region data associated with gene IL4R ...
[2022-08-18 14:08:47.542644 EDT] Extracting fasta sequence associated with chromosome: 16, tss: 27313974, strand: 1, gene: IL4R ...
[2022-08-18 14:08:47.542673 EDT] Processing mapped ambiguity code TGC ...
[2022-08-18 14:08:47.542684 EDT] Processing mapped ambiguity code TAC ...
[2022-08-18 14:08:47.542692 EDT] Processing mapped ambiguity code AGC ...
[2022-08-18 14:08:47.542700 EDT] Processing mapped ambiguity code AAC ...
[2022-08-18 14:08:47.542710 EDT] Processing region data associated with gene AFF3 ...
[2022-08-18 14:08:47.542718 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and AFF3 strand orientation is -1 ...
[2022-08-18 14:08:47.542737 EDT] Processing region data associated with gene CXCR4 ...
[2022-08-18 14:08:47.542747 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and CXCR4 strand orientation is -1 ...
[2022-08-18 14:08:47.542764 EDT] Processing region data associated with gene GNA13 ...
[2022-08-18 14:08:47.542804 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and GNA13 strand orientation is -1 ...
[2022-08-18 14:08:47.542823 EDT] Processing region data associated with gene IGLL5 ...
[2022-08-18 14:08:47.542833 EDT] Extracting fasta sequence associated with chromosome: 22, tss: 22887816, strand: 1, gene: IGLL5 ...
[2022-08-18 14:08:47.542844 EDT] Processing mapped ambiguity code TGC ...
[2022-08-18 14:08:47.542852 EDT] Processing mapped ambiguity code TAC ...
[2022-08-18 14:08:47.542861 EDT] Processing mapped ambiguity code AGC ...
[2022-08-18 14:08:47.542896 EDT] Processing mapped ambiguity code AAC ...
[2022-08-18 14:08:47.542923 EDT] Processing region data associated with gene BCL7A ...
[2022-08-18 14:08:47.542935 EDT] Extracting fasta sequence associated with chromosome: 12, tss: 122019422, strand: 1, gene: BCL7A ...
[2022-08-18 14:08:47.542946 EDT] Processing mapped ambiguity code TGC ...
[2022-08-18 14:08:47.542967 EDT] Processing mapped ambiguity code TAC ...
[2022-08-18 14:08:47.542996 EDT] Processing mapped ambiguity code AGC ...
[2022-08-18 14:08:47.543005 EDT] Processing mapped ambiguity code AAC ...
[2022-08-18 14:08:47.543025 EDT] Processing region data associated with gene CD83 ...
[2022-08-18 14:08:47.543034 EDT] Extracting fasta sequence associated with chromosome: 6, tss: 14117256, strand: 1, gene: CD83 ...
[2022-08-18 14:08:47.543059 EDT] Processing mapped ambiguity code TGC ...
[2022-08-18 14:08:47.543070 EDT] Processing mapped ambiguity code TAC ...
[2022-08-18 14:08:47.543078 EDT] Processing mapped ambiguity code AGC ...
[2022-08-18 14:08:47.543086 EDT] Processing mapped ambiguity code AAC ...
[2022-08-18 14:08:47.543097 EDT] Processing region data associated with gene H2AC16 ...
[2022-08-18 14:08:47.543105 EDT] Extracting fasta sequence associated with chromosome: 6, tss: 27865317, strand: 1, gene: H2AC16 ...
[2022-08-18 14:08:47.543115 EDT] Processing mapped ambiguity code TGC ...
[2022-08-18 14:08:47.543123 EDT] Processing mapped ambiguity code TAC ...
[2022-08-18 14:08:47.543133 EDT] Processing mapped ambiguity code AGC ...
[2022-08-18 14:08:47.543142 EDT] Processing mapped ambiguity code AAC ...
[2022-08-18 14:08:47.543152 EDT] Processing region data associated with gene TP53 ...
[2022-08-18 14:08:47.543160 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and TP53 strand orientation is -1 ...
[2022-08-18 14:08:47.553117 EDT] Processing region data associated with gene SOCS1 ...
[2022-08-18 14:08:47.553132 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and SOCS1 strand orientation is -1 ...
[2022-08-18 14:08:47.553179 EDT] Processing region data associated with gene BTG1 ...
[2022-08-18 14:08:47.553192 EDT] Extracting fasta sequence associated with chromosome: 12, tss: 92145846, strand: -1, gene: BTG1 ...
[2022-08-18 14:08:47.553251 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.553274 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.553290 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.553305 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.553319 EDT] Processing region data associated with gene MYH10 ...
[2022-08-18 14:08:47.553328 EDT] Extracting fasta sequence associated with chromosome: 17, tss: 8630761, strand: -1, gene: MYH10 ...
[2022-08-18 14:08:47.553338 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.553346 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.557015 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.557030 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.557058 EDT] Processing region data associated with gene ARID5B ...
[2022-08-18 14:08:47.557073 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and ARID5B strand orientation is 1 ...
[2022-08-18 14:08:47.557090 EDT] Processing region data associated with gene EBF1 ...
[2022-08-18 14:08:47.557111 EDT] Extracting fasta sequence associated with chromosome: 5, tss: 159099916, strand: -1, gene: EBF1 ...
[2022-08-18 14:08:47.557122 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.557140 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.557149 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.580097 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.580135 EDT] Processing region data associated with gene H1-5 ...
[2022-08-18 14:08:47.580146 EDT] Extracting fasta sequence associated with chromosome: 6, tss: 27867588, strand: -1, gene: H1-5 ...
[2022-08-18 14:08:47.580159 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.580184 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.580201 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.580211 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.580222 EDT] Processing region data associated with gene H1-2 ...
[2022-08-18 14:08:47.580231 EDT] Extracting fasta sequence associated with chromosome: 6, tss: 26056470, strand: -1, gene: H1-2 ...
[2022-08-18 14:08:47.580243 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.580252 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.580261 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.580270 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.580279 EDT] Processing region data associated with gene SGK1 ...
[2022-08-18 14:08:47.580289 EDT] Extracting fasta sequence associated with chromosome: 6, tss: 134318112, strand: -1, gene: SGK1 ...
[2022-08-18 14:08:47.580310 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.580321 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.580330 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.580340 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.580348 EDT] Processing region data associated with gene IL4R ...
[2022-08-18 14:08:47.580359 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and IL4R strand orientation is 1 ...
[2022-08-18 14:08:47.580379 EDT] Processing region data associated with gene AFF3 ...
[2022-08-18 14:08:47.580388 EDT] Extracting fasta sequence associated with chromosome: 2, tss: 100106128, strand: -1, gene: AFF3 ...
[2022-08-18 14:08:47.580426 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.580439 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.580449 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.580459 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.580468 EDT] Processing region data associated with gene CXCR4 ...
[2022-08-18 14:08:47.580476 EDT] Extracting fasta sequence associated with chromosome: 2, tss: 136116243, strand: -1, gene: CXCR4 ...
[2022-08-18 14:08:47.580486 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.580494 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.580501 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.580509 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.580516 EDT] Processing region data associated with gene GNA13 ...
[2022-08-18 14:08:47.580525 EDT] Extracting fasta sequence associated with chromosome: 17, tss: 65056740, strand: -1, gene: GNA13 ...
[2022-08-18 14:08:47.580535 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.580543 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.580589 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.580600 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.580609 EDT] Processing region data associated with gene IGLL5 ...
[2022-08-18 14:08:47.580617 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and IGLL5 strand orientation is 1 ...
[2022-08-18 14:08:47.580633 EDT] Processing region data associated with gene BCL7A ...
[2022-08-18 14:08:47.580641 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and BCL7A strand orientation is 1 ...
[2022-08-18 14:08:47.580657 EDT] Processing region data associated with gene CD83 ...
[2022-08-18 14:08:47.580665 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and CD83 strand orientation is 1 ...
[2022-08-18 14:08:47.580695 EDT] Processing region data associated with gene H2AC16 ...
[2022-08-18 14:08:47.580708 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and H2AC16 strand orientation is 1 ...
[2022-08-18 14:08:47.580723 EDT] Processing region data associated with gene TP53 ...
[2022-08-18 14:08:47.580733 EDT] Extracting fasta sequence associated with chromosome: 17, tss: 7687491, strand: -1, gene: TP53 ...
[2022-08-18 14:08:47.580743 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.580751 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.580778 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.580790 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.580799 EDT] Processing region data associated with gene SOCS1 ...
[2022-08-18 14:08:47.580858 EDT] Extracting fasta sequence associated with chromosome: 16, tss: 11256200, strand: -1, gene: SOCS1 ...
[2022-08-18 14:08:47.580873 EDT] Processing mapped ambiguity code GCA ...
[2022-08-18 14:08:47.580881 EDT] Processing mapped ambiguity code GTT ...
[2022-08-18 14:08:47.580889 EDT] Processing mapped ambiguity code GTA ...
[2022-08-18 14:08:47.580909 EDT] Processing mapped ambiguity code GCT ...
[2022-08-18 14:08:47.580919 EDT] Preparing ambiguity code strings final analysis ...
[2022-08-18 14:08:47.581058 EDT] Preparing variants final analysis ...
[2022-08-18 14:08:47.663842 EDT] Preparing to print output CSV files ...
[2022-08-18 14:08:47.663904 EDT] Printing output CSV files ...
[2022-08-18 14:08:47.674745 EDT] Shutting down Fasta Region Inspector v0.1.0.0 ...
```

## Dependency Visualization via Stack Dot
The following is a visualization of the dependency graph for **FRI**:


![](dependencies.png)
