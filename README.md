# Fasta-Region-Inspector

## Introduction
Fasta-Region-Inspector (**FRI**) is a bioinformatics tool for analyzing somatic hypermutation (SHM).

## Purpose
Analyzing cancer cohorts for SHM is a interesting cancer-genomics related research question, and one that is computationally intensive.  SHM occurs when immunoglobulin genes accumulate apparently random point mutations within productively rearranged V, D, and J segments, as defined by [Elaine S. Jaffe MD, in Hematopathology, 2017](https://www.sciencedirect.com/topics/immunology-and-microbiology/somatic-hypermutation).  More specifically, [Blood. 2009 Apr 16; 113(16): 3706–3715](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2670789/) states that the hallmark of SHM is the increased percentage of mutations in hypermutable RGYW (underlined residue is the most frequently mutated) motifs.  Correctly identifying SHM in large cancer genomics datasets can be a complex and tedious computational problem. 

This bioinformatics tool gives researchers the ability to answer the following questions:
- Which variants are within 2 kb of the transcription start site (TSS) of the corresponding gene?
- Where do user-defined mapped ambiguity strings lie within 2 kb of the TSS?
- Amalgamating the previous 2 questions, which variant(s) are found within a mapped ambiguity string, that also lie within 2 kb of the TSS?

This tool aims to answer common SHM variant-level questions in a software package that provides:
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
  - After all of the above changes (moving to stack or the new YAML input file doesn't affect runtime performance, but the rest of the points do), there is ~160X speedup in runtime (benchmarks up-coming) in the examples run locally so far (comparable tests are showing ~8 hrs down to ~3 minutes).
  - In general, the more ambiguity codes and variants there are to be analyzed by **FRI**, the larger the runtime disparity between the old and new version.

### Conclusion
Don't use the old version!

## Configuration YAML

**FRI** utilizes a configuration YAML file to provide the necessary components for a successful run.

The following keys are **required**:
- ```Fasta``` -> The filepath to the FASTA file. (String)
- ```Variants``` -> A [compact nested mapping](https://yaml.org/spec/1.2.2/#chapter-2-language-overview) (see [below](https://github.com/Matthew-Mosior/Fasta-Region-Inspector/blob/main/README.md#variant-type) for more information).
- ```Ambiguity_Codes``` -> An array of ambiguity codes to search against. (String)
- ```Output_Directory``` -> The filepath to the output directory (must already exist). (String)
- ```Keep_BioMart``` -> Whether or not to keep the BioMart file from the wget system call. (Boolean)
- ```Ignore_Strandedness``` -> Whether or not to ignore strandness of the respective gene and search the TSS in both directions. (Boolean)

The following keys are **optional**:
- ```TSS_Window_Size``` -> The TSS window size to search across. (String)

## Variant Type
The ```Variant``` [compact nested mapping](https://yaml.org/spec/1.2.2/#chapter-2-language-overview) represents the required information surrounding a variant-of-interest to be examined by **FRI**.

The following keys are **required**:
- ```Sample``` -> The associated sample identifier for the respective variant (String).
- ```HGNC_Symbol``` -> The HGNC symbol for the respective variant (String).
- ```Chromosome``` -> The chromosome the respective variant lies within (String).
- ```Start_Position``` -> The start position for the respective variant (String).
- ```End_Position``` -> The end position for the respective variant (String).
- ```Reference_Allele``` -> The reference allele for the respective variant (String).
- ```Alternate_Allele``` -> The alternate allele for the respective variant (String).
- ```ENST``` -> The ENST for the respective variant (String).

### Note
The above data can typically be easily grepped/programmed for from the output of a bioinformatics pipeline run, such as the [alignment_exome.cwl](https://github.com/genome/analysis-workflows/blob/master/definitions/pipelines/alignment_exome.cwl).

You can also start from a [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)-annotated variant file, which can typically be created once you have run an aligner/variant caller (see above).

## Building the project
This software was developed on a M1 mac.  It has not been tested on other operating systems or chipsets, but should run perfectly fine on any compatible platform.

If you are on a M1 mac, please build the project using the following:
```
% stack build --arch aarch64
```

If you are **NOT** on a M1 mac, please build the project using the following:
```
% stack build
```

## Example Usage
**FRI** is easy to use, as it only requires a single command-line positional argument, the configuration YAML.

If you are on a M1 mac, please run the project using the following:
```
% stack exec --arch aarch64 fasta-region-inspector-exe /path/to/configuration.yaml
```

If you are **NOT** a M1 mac, please run the project using the following:
```
% stack exec fasta-region-inspector-exe /path/to/configuration.yaml
```


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

[2022-08-15 17:03:08.951069 EDT] Starting up Fasta Region Inspector v0.1.0.0 ...
[2022-08-15 17:03:08.952688 EDT] Query BioMart for regions data ...
[2022-08-15 17:03:08.952745 EDT] Generating BioMart compatible XML ...
[2022-08-15 17:03:08.952767 EDT] Querying and downloading region data from BioMart via system process call ...
--2022-08-15 17:03:09--  http://www.ensembl.org/biomart/martservice?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20count=%22%22%20datasetConfigVersion=%220.6%22%20formatter=%22CSV%22%20header=%220%22%20uniqueRows=%220%22%20virtualSchemaName=%22default%22%3E%3CDataset%20interface=%22default%22%20name=%22hsapiens_gene_ensembl%22%3E%3CFilter%20name=%22ensembl_transcript_id%22%20value=%22ENST00000409579,ENST00000279873,ENST00000538010,ENST00000256015,ENST00000612003,ENST00000409817,ENST00000313708,ENST00000439174,ENST00000331442,ENST00000613174,ENST00000343677,ENST00000526893,ENST00000395762,ENST00000269243,ENST00000367858,ENST00000644787,ENST00000635293%22/%3E%3CAttribute%20name=%22chromosome_name%22/%3E%3CAttribute%20name=%22transcription_start_site%22/%3E%3CAttribute%20name=%22strand%22/%3E%3CAttribute%20name=%22external_gene_name%22/%3E%3C/Dataset%3E%3C/Query%3E
Resolving www.ensembl.org (www.ensembl.org)... 193.62.193.83
Connecting to www.ensembl.org (www.ensembl.org)|193.62.193.83|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: unspecified [text/plain]
Saving to: ‘/Users/matthewmosior/Software/local/haskell/Fasta-Region-Inspector-Updated/Compact-Regions-Implementation-New/examples/test-output/biomartresult.txt’

/Users/matthewmosior/Software/local/haskell     [ <=>                                                                                     ]     339  --.-KB/s    in 0s

2022-08-15 17:03:09 (7.18 MB/s) - ‘/Users/matthewmosior/Software/local/haskell/Fasta-Region-Inspector-Updated/Compact-Regions-Implementation-New/examples/test-output/biomartresult.txt’ saved [339]

[2022-08-15 17:03:09.560481 EDT] Successfully queried and downloaded region data from BioMart via system process call ...
[2022-08-15 17:03:09.560623 EDT] Beginning to process downloaded region data ...
[2022-08-15 17:03:09.560864 EDT] Determining whether each variant is within the its respective genes TSS ...
[2022-08-15 17:03:09.560910 EDT] Calculating the reverse complement of each user defined ambiguity code ...
[2022-08-15 17:03:09.560956 EDT] Creating list of tuples to define directionality of each forward strand ambiguity code ...
[2022-08-15 17:03:09.560997 EDT] Creating list of tuples to define directionality of each reverse strand ambiguity code ...
[2022-08-15 17:03:09.561031 EDT] Generating all possible ambiguity codes strings using SMT solver ...
[2022-08-15 17:03:09.657130 EDT] Preparing ambiguity code strings to determine whether each lies within its respective TSS ...
[2022-08-15 17:03:09.657179 EDT] Reading fasta file into strict ByteString ...
[2022-08-15 17:03:12.311377 EDT] Parsing fasta file into [Sequence] ...
[2022-08-15 17:04:06.339396 EDT] Extract every sequences name and characters from cfasta ...
[2022-08-15 17:04:06.340851 EDT] Putting namesandcharacters into compact region ...
[2022-08-15 17:06:00.422136 EDT] Determing whether each ambiguity code string lies within its respective TSS ...
[2022-08-15 17:06:00.441892 EDT] Processing region data associated with gene BTG1 ...
[2022-08-15 17:06:00.441979 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and BTG1 strand orientation is -1 ...
[2022-08-15 17:06:00.442028 EDT] Processing region data associated with gene MYH10 ...
[2022-08-15 17:06:00.442044 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and MYH10 strand orientation is -1 ...
[2022-08-15 17:06:00.442068 EDT] Processing region data associated with gene ARID5B ...
[2022-08-15 17:06:00.442089 EDT] Extracting fasta sequence associated with chromosome: 10, tss: 61901699, strand: 1, gene: ARID5B ...
[2022-08-15 17:06:00.442109 EDT] Processing mapped ambiguity code TGC ...
[2022-08-15 17:06:00.447462 EDT] Processing mapped ambiguity code TAC ...
[2022-08-15 17:06:00.447592 EDT] Processing mapped ambiguity code AGC ...
[2022-08-15 17:06:00.447644 EDT] Processing mapped ambiguity code AAC ...
[2022-08-15 17:06:00.447681 EDT] Processing region data associated with gene EBF1 ...
[2022-08-15 17:06:00.447702 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and EBF1 strand orientation is -1 ...
[2022-08-15 17:06:00.458987 EDT] Processing region data associated with gene H1-5 ...
[2022-08-15 17:06:00.459031 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and H1-5 strand orientation is -1 ...
[2022-08-15 17:06:00.464630 EDT] Processing region data associated with gene H1-2 ...
[2022-08-15 17:06:00.464673 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and H1-2 strand orientation is -1 ...
[2022-08-15 17:06:00.464740 EDT] Processing region data associated with gene SGK1 ...
[2022-08-15 17:06:00.467068 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and SGK1 strand orientation is -1 ...
[2022-08-15 17:06:00.467250 EDT] Processing region data associated with gene IL4R ...
[2022-08-15 17:06:00.470870 EDT] Extracting fasta sequence associated with chromosome: 16, tss: 27313974, strand: 1, gene: IL4R ...
[2022-08-15 17:06:00.470965 EDT] Processing mapped ambiguity code TGC ...
[2022-08-15 17:06:00.471037 EDT] Processing mapped ambiguity code TAC ...
[2022-08-15 17:06:00.471125 EDT] Processing mapped ambiguity code AGC ...
[2022-08-15 17:06:00.471153 EDT] Processing mapped ambiguity code AAC ...
[2022-08-15 17:06:00.471173 EDT] Processing region data associated with gene AFF3 ...
[2022-08-15 17:06:00.471184 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and AFF3 strand orientation is -1 ...
[2022-08-15 17:06:00.471216 EDT] Processing region data associated with gene CXCR4 ...
[2022-08-15 17:06:00.480915 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and CXCR4 strand orientation is -1 ...
[2022-08-15 17:06:00.480985 EDT] Processing region data associated with gene GNA13 ...
[2022-08-15 17:06:00.481002 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and GNA13 strand orientation is -1 ...
[2022-08-15 17:06:00.501160 EDT] Processing region data associated with gene IGLL5 ...
[2022-08-15 17:06:00.501208 EDT] Extracting fasta sequence associated with chromosome: 22, tss: 22887816, strand: 1, gene: IGLL5 ...
[2022-08-15 17:06:00.501224 EDT] Processing mapped ambiguity code TGC ...
[2022-08-15 17:06:00.509267 EDT] Processing mapped ambiguity code TAC ...
[2022-08-15 17:06:00.509313 EDT] Processing mapped ambiguity code AGC ...
[2022-08-15 17:06:00.509324 EDT] Processing mapped ambiguity code AAC ...
[2022-08-15 17:06:00.509347 EDT] Processing region data associated with gene BCL7A ...
[2022-08-15 17:06:00.509362 EDT] Extracting fasta sequence associated with chromosome: 12, tss: 122019422, strand: 1, gene: BCL7A ...
[2022-08-15 17:06:00.509377 EDT] Processing mapped ambiguity code TGC ...
[2022-08-15 17:06:00.509388 EDT] Processing mapped ambiguity code TAC ...
[2022-08-15 17:06:00.509396 EDT] Processing mapped ambiguity code AGC ...
[2022-08-15 17:06:00.509406 EDT] Processing mapped ambiguity code AAC ...
[2022-08-15 17:06:00.509416 EDT] Processing region data associated with gene CD83 ...
[2022-08-15 17:06:00.509425 EDT] Extracting fasta sequence associated with chromosome: 6, tss: 14117256, strand: 1, gene: CD83 ...
[2022-08-15 17:06:00.509436 EDT] Processing mapped ambiguity code TGC ...
[2022-08-15 17:06:00.509452 EDT] Processing mapped ambiguity code TAC ...
[2022-08-15 17:06:00.509460 EDT] Processing mapped ambiguity code AGC ...
[2022-08-15 17:06:00.509471 EDT] Processing mapped ambiguity code AAC ...
[2022-08-15 17:06:00.509482 EDT] Processing region data associated with gene H2AC16 ...
[2022-08-15 17:06:00.509551 EDT] Extracting fasta sequence associated with chromosome: 6, tss: 27865317, strand: 1, gene: H2AC16 ...
[2022-08-15 17:06:00.509573 EDT] Processing mapped ambiguity code TGC ...
[2022-08-15 17:06:00.509598 EDT] Processing mapped ambiguity code TAC ...
[2022-08-15 17:06:00.509611 EDT] Processing mapped ambiguity code AGC ...
[2022-08-15 17:06:00.509622 EDT] Processing mapped ambiguity code AAC ...
[2022-08-15 17:06:00.509639 EDT] Processing region data associated with gene TP53 ...
[2022-08-15 17:06:00.509653 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and TP53 strand orientation is -1 ...
[2022-08-15 17:06:00.509676 EDT] Processing region data associated with gene SOCS1 ...
[2022-08-15 17:06:00.509686 EDT] Could not process region data associated with current ambiguity code WRC:
                                 WRC strand orientation is 1 and SOCS1 strand orientation is -1 ...
[2022-08-15 17:06:00.509753 EDT] Processing region data associated with gene BTG1 ...
[2022-08-15 17:06:00.509772 EDT] Extracting fasta sequence associated with chromosome: 12, tss: 92145846, strand: -1, gene: BTG1 ...
[2022-08-15 17:06:00.509787 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.509812 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.509850 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.509956 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.509991 EDT] Processing region data associated with gene MYH10 ...
[2022-08-15 17:06:00.510109 EDT] Extracting fasta sequence associated with chromosome: 17, tss: 8630761, strand: -1, gene: MYH10 ...
[2022-08-15 17:06:00.510191 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.510212 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.510225 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.510233 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.510242 EDT] Processing region data associated with gene ARID5B ...
[2022-08-15 17:06:00.510250 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and ARID5B strand orientation is 1 ...
[2022-08-15 17:06:00.510271 EDT] Processing region data associated with gene EBF1 ...
[2022-08-15 17:06:00.510282 EDT] Extracting fasta sequence associated with chromosome: 5, tss: 159099916, strand: -1, gene: EBF1 ...
[2022-08-15 17:06:00.510295 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.510303 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.510311 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.510320 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.510332 EDT] Processing region data associated with gene H1-5 ...
[2022-08-15 17:06:00.510371 EDT] Extracting fasta sequence associated with chromosome: 6, tss: 27867588, strand: -1, gene: H1-5 ...
[2022-08-15 17:06:00.510438 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.510459 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.510474 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.510484 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.510494 EDT] Processing region data associated with gene H1-2 ...
[2022-08-15 17:06:00.510617 EDT] Extracting fasta sequence associated with chromosome: 6, tss: 26056470, strand: -1, gene: H1-2 ...
[2022-08-15 17:06:00.510638 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.510650 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.510659 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.510667 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.510678 EDT] Processing region data associated with gene SGK1 ...
[2022-08-15 17:06:00.510688 EDT] Extracting fasta sequence associated with chromosome: 6, tss: 134318112, strand: -1, gene: SGK1 ...
[2022-08-15 17:06:00.510711 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.510763 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.510773 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.511783 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.511797 EDT] Processing region data associated with gene IL4R ...
[2022-08-15 17:06:00.511807 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and IL4R strand orientation is 1 ...
[2022-08-15 17:06:00.511895 EDT] Processing region data associated with gene AFF3 ...
[2022-08-15 17:06:00.511939 EDT] Extracting fasta sequence associated with chromosome: 2, tss: 100106128, strand: -1, gene: AFF3 ...
[2022-08-15 17:06:00.511956 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.511974 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.512059 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.512078 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.512090 EDT] Processing region data associated with gene CXCR4 ...
[2022-08-15 17:06:00.512104 EDT] Extracting fasta sequence associated with chromosome: 2, tss: 136116243, strand: -1, gene: CXCR4 ...
[2022-08-15 17:06:00.516016 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.516083 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.516098 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.516109 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.516121 EDT] Processing region data associated with gene GNA13 ...
[2022-08-15 17:06:00.516155 EDT] Extracting fasta sequence associated with chromosome: 17, tss: 65056740, strand: -1, gene: GNA13 ...
[2022-08-15 17:06:00.516170 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.516181 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.516191 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.516288 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.516347 EDT] Processing region data associated with gene IGLL5 ...
[2022-08-15 17:06:00.516361 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and IGLL5 strand orientation is 1 ...
[2022-08-15 17:06:00.516381 EDT] Processing region data associated with gene BCL7A ...
[2022-08-15 17:06:00.516398 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and BCL7A strand orientation is 1 ...
[2022-08-15 17:06:00.516445 EDT] Processing region data associated with gene CD83 ...
[2022-08-15 17:06:00.516459 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and CD83 strand orientation is 1 ...
[2022-08-15 17:06:00.516912 EDT] Processing region data associated with gene H2AC16 ...
[2022-08-15 17:06:00.516941 EDT] Could not process region data associated with current ambiguity code GYW:
                                 GYW strand orientation is -1 and H2AC16 strand orientation is 1 ...
[2022-08-15 17:06:00.517000 EDT] Processing region data associated with gene TP53 ...
[2022-08-15 17:06:00.517020 EDT] Extracting fasta sequence associated with chromosome: 17, tss: 7687491, strand: -1, gene: TP53 ...
[2022-08-15 17:06:00.517038 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.517060 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.517089 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.517098 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.517122 EDT] Processing region data associated with gene SOCS1 ...
[2022-08-15 17:06:00.517525 EDT] Extracting fasta sequence associated with chromosome: 16, tss: 11256200, strand: -1, gene: SOCS1 ...
[2022-08-15 17:06:00.517543 EDT] Processing mapped ambiguity code GCA ...
[2022-08-15 17:06:00.517554 EDT] Processing mapped ambiguity code GTT ...
[2022-08-15 17:06:00.517612 EDT] Processing mapped ambiguity code GTA ...
[2022-08-15 17:06:00.517626 EDT] Processing mapped ambiguity code GCT ...
[2022-08-15 17:06:00.517638 EDT] Preparing ambiguity code strings final analysis ...
[2022-08-15 17:06:00.518505 EDT] Preparing variants final analysis ...
[2022-08-15 17:06:00.527625 EDT] Preparing to print output CSV files ...
[2022-08-15 17:06:00.527921 EDT] Printing output CSV files ...
[2022-08-15 17:06:00.575115 EDT] Shutting down Fasta Region Inspector v0.1.0.0 ...
```
