# fasta-region-inspector

## Introduction
fasta-region-inspector (**FRI**) is a bioinformatics tool for analyzing somatic hypermutation (SHM).

## Purpose
Analyzing cancer cohorts for SHM is a interesting cancer-genomics related research question, and one that is computationally intensive.  SHM occurs when immunoglobulin genes accumulate apparently random point mutations within productively rearranged V, D, and J segments, as defined by [Elaine S. Jaffe MD, in Hematopathology, 2017](https://www.sciencedirect.com/topics/immunology-and-microbiology/somatic-hypermutation).  More specifically, [Blood. 2009 Apr 16; 113(16): 3706–3715](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2670789/) states that the hallmark of SHM is the increased percentage of mutations in hypermutable RGYW motifs.  Correctly identifying SHM in large cancer genomics datasets can be a complex and tedious computational problem. 

This bioinformatics tool gives researchers the ability to answer the following questions:
- Which variants are within 2 kb of the transcription start site (TSS) of the corresponding gene?
- Where do user-defined mapped ambiguity strings lie within 2 kb of the TSS?
- Amalgamating the previous 2 questions, which variant(s) are found within a mapped ambiguity string, that also lie within 2 kb of the TSS?

## Improvements from Previous Implementation

### Overview
This tool aims to answer common SHM variant-level questions in a software package that provides:
- Excellent runtime performance in a robust, functional implementation.
- Minimized memory usage.
- A simple, YAML input file format.
- Clean, informative stdout logging.

### Old vs. New
This version of **FRI** is vastly improved upon in numerous ways, including:
- Optimized file IO via `linear types`.
  - In the previous version, the entire input FASTA file was read into memory, irrespective of the actual needs for the specific analysis.
  - Using `linear types`, system resource usage is provably safe (file handles must be used **once**, and only once).
  - This allows for an incredibly efficient, highly tuned algorithm to seek as close to the TSS sequence as possible, without wasting memory or time.
  - This new and vastly improved method of file IO is a large contributor to the large runtime disparity between the old and new version (**~8X** speedup).
- Much more efficient use of multi-core CPU architecture via multi-threading.
  - The previous version lacked any sense of multi-threading.
  - The new version utilizes scope-delimited thread lifetimes in order to concurrently perform the most computationally expensive functions in parallel.
- A new effect system.
  - The previous version lacked any sort of effect system, mtl or otherwise.
  - The new version introduces effects via the `effectful` ecosystem, which is a series of packages optimized for performance, usability, and clarity.
- More in-depth logging.
  - The previous version handled logging as a simple call to the clock, and a statement.
  - The new version handles logging in a more comprehensive fashion, incorporating thread information (due to the newly introduced multi-threaded computational model, see above), as well as the originating function which contains the logging statement, the level of logging, and finally the current time via the clock.
- General code clean-up
  - The plethora of default-extensions in the `package.yaml` was removed in favor of a more direct and streamline approach for language extension handling, wherein each individual module contains only the extensions it needs.
  - Many functions/functionality were re-organized into a wider and thinner module structure (slightly more modules, and less code per module) to help with general organization.

## Configuration YAML

**FRI** utilizes a configuration YAML file to provide the necessary components for a successful run.

The following keys are **required**:
- ```Fasta``` -> The filepath to the FASTA file. (String)
- ```Fasta_Index``` -> The filepath to the FASTA Index (fai) file. (String)
- ```Variants``` -> A [compact nested mapping](https://yaml.org/spec/1.2.2/#chapter-2-language-overview) (see [below](https://github.com/Matthew-Mosior/Fasta-Region-Inspector/blob/main/README.md#variant-type) for more information).
- ```Ambiguity_Codes``` -> An array of ambiguity codes to search against. (String)
- ```Output_Directory``` -> The filepath to the output directory (must already exist). (String)
- ```Keep_BioMart``` -> Whether or not to keep the BioMart file from the http request. (Boolean)
- ```Ignore_Strandedness``` -> Whether or not to ignore strandness of the respective gene and search the TSS in both directions. (Boolean)
- ```Write_Ambiguity_Codes``` -> Whether or not to write CSV file containing all locations of user-defined mapped ambiguity code strings in TSS window. (Boolean)
- ```Max_Number_of_Concurrent_Threads``` -> The max number of threads to be able to use at any one time. (Int)
  - You can use up to `999` threads.

The following keys are **optional**:
- ```TSS_Window_Size``` -> The TSS window size to search across. (String)
  - The default is 2000 bp.

## Variant Type
The ```Variants``` [compact nested mapping](https://yaml.org/spec/1.2.2/#chapter-2-language-overview) key represents the required information surrounding a variant-of-interest to be examined by **FRI**.

The following keys are **required**:
- ```Sample``` -> The associated sample identifier for the respective variant (String).
- ```HGNC_Symbol``` -> The HGNC symbol for the respective variant (String).
- ```Chromosome``` -> The chromosome the respective variant lies within (String).
  - You should include all characters up to the first tab character in the description line.
  - For example, given this sequence description line, ```>2 dna:chromosome chromosome:GRCh38:2:1:242193529:1 REF```, you would put ```>2```.
- ```Start_Position``` -> The start position for the respective variant (String).
- ```End_Position``` -> The end position for the respective variant (String).
- ```Reference_Allele``` -> The reference allele for the respective variant (String).
- ```Alternate_Allele``` -> The alternate allele for the respective variant (String).
- ```ENST``` -> The ENST for the respective variant (String).

### Notes
The above data can typically be easily grepped/programmed for from the output of a bioinformatics pipeline run, such as the [alignment_exome.cwl](https://github.com/genome/analysis-workflows/blob/master/definitions/pipelines/alignment_exome.cwl).

You can also start from a [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html)-annotated variant file, which can typically be created once you have run an aligner/variant caller (see above).  Starting from a VEP-annotated file will allow you to filter the variants on gene so that you can better target your analysis with this tool.

## Dependency - Theorem Solver
The program uses the [sbv](https://hackage.haskell.org/package/sbv) package, which interfaces with theorem solver(s).

This package interfaces with one of following supported theorem solvers:
- [ABC](http://www.eecs.berkeley.edu/~alanmi/abc/) from University of Berkeley.
- [Boolector](http://fmv.jku.at/boolector/) from Johannes Kepler University.
- [Bitwuzla](http://bitwuzla.github.io/) from Stanford University.
- [CVC4](http://cvc4.github.io/) and [CVC5](http://cvc5.github.io/) from Stanford University and the University of Iowa.
- [DReal](http://dreal.github.com/) from CMU.
- [MathSAT](http://mathsat.fbk.eu/) from FBK and DISI-University of Trento.
- [Yices](http://yices.csl.sri.com/) from SRI.
- [Z3](http://github.com/Z3Prover/z3/) from Microsoft.

This software has been tested with both Z3 and Yices, although the other theorem solvers listed above should work just as well.

You will need to install and set-up one of the above theorem solvers in order for **FRI** to calculate the mapped ambiguity code string for the given ambiguity code(s).

Please see http://leventerkok.github.io/sbv/ to get more information on how theorem solvers work.

## Dependency - LLVM
**FRI** is set up to compile via [LLVM](https://llvm.org/) instead of using the native code generator.

This will generally take slightly longer than the native code generator to compile. Produced code is generally the same speed or faster than the other two code generators. Compiling via LLVM requires LLVM’s opt and llc executables to be in PATH (taken from GHC documentation).

This means you need the LLVM compiler to compile this code.

If you are on a mac, you can easily install LLVM via brew:

```brew install llvm```

If you are on windows, you can follow [this](https://llvm.org/docs/GettingStartedVS.html) guide to get started with LLVM.

If you are linux/unix, please download and install via your distribution's package manager.

The `package.yaml` file is where this option is set via `ghc-options` using the `-fllvm` flag.

## Building the project
This software was developed on an M1 mac.

You will need:
- The GHC compiler
  - I recommend downloading this via [GHCup](https://www.haskell.org/ghcup/).
- The Haskell Stack Tool (https://docs.haskellstack.org/en/stable/install_and_upgrade/).

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
% stack exec --arch aarch64 fasta-region-inspector /path/to/configuration.yaml
```

If you are **NOT** an M1 mac, please run the project using the following:
```
% stack exec fasta-region-inspector /path/to/configuration.yaml
```

## Output
**By default**, **FRI** outputs a single CSV file,namely *variants_in_ambiguity_codes.csv*, containing all of the information needed to answer the questions posed at the top of this README.

Each user-defined variant will be included in this CSV file, with the following fields:
- ```Variant``` -> A description of the user-defined variant.
- ```Region``` -> A description of the region related metadata pulled from BioMart for the user-defined variant.
- ```Variant_Within_Region``` -> A boolean (Y/N) of whether or not the user-defined variant lies within 2 kb of the TSS.
- ```Ambiguity_Code``` -> The user-defined ambiguity code used for the stringsearch.
  - This can be the reverse complement of the user-defined ambiguity code, depending on the strand directionality pulled from BioMart for the variant.
- ```Mapped_Nucleotide_String``` -> The [Satisfiability Modulo Theories](https://hackage.haskell.org/package/sbv) calculated mapped ambuguity code string for the given ambiguity code.
- ```Ambiguity_Code_String_Locations_Within_TSS``` -> The starting position of the first base pair for the mapped nucleotide string.
  - Again, the strand directionality can change which way the mapped nucleotide string is interpreted (see below).

### Notes

```N/A``` values in the ```Ambiguity_Code```, ```Mapped_Nucleotide_String```, or ```Ambiguity_Code_String_Locations_Within_TSS``` fields indicate that the respective variant did not lie within a mapped nucleotide string within 2 kb of the TSS.

Depending on how you set the ```Write_Ambiguity_Codes``` key in the configuration YAML, you can also get an additional CSV file containing all locations of user-defined mapped ambiguity code strings in the TSS window.

## Visualizing the Output
Interpretation of the *variants_in_ambiguity_codes.csv* file is a simple process.

The following example will hopefully help with understanding the output of **FRI**.

### Note

If you wish to follow along with your own data, please download [IGV](https://software.broadinstitute.org/software/igv/).

IGV is a helpful resource for visualizing your results (all of the following screenshots in this section are from IGV). 

### Examine the Mapped Nucleotide String<sup>**Assumes a working knowledge of IGV**</sup>

The following assumes you load the FASTA file used for the **FRI** run.

|Variant|Region|Variant_Within_Region|Ambiguity_Code|Mapped_Nucleotide_String|Ambiguity_Code_String_Locations_Within_TSS|
|-------|------|---------------------|--------------|------------------------|------------------------------------------|
|SMP1:H1-5:chr6:27867219:27867219:G:A:ENST00000331442|6:27867588:&#x200B;-1:H1-5|Y|GYW|GCT|27867220|

Let's take a look at what IGV shows for this particular variant:
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
                                         ____           __                              _                   _                            __            
                                        / __/___ ______/ /_____ _      ________  ____ _(_)___  ____        (_)___  _________  ___  _____/ /_____  _____
                                       / /_/ __ `/ ___/ __/ __ `/_____/ ___/ _ \/ __ `/ / __ \/ __ \______/ / __ \/ ___/ __ \/ _ \/ ___/ __/ __ \/ ___/
                                      / __/ /_/ (__  ) /_/ /_/ /_____/ /  /  __/ /_/ / / /_/ / / / /_____/ / / / (__  ) /_/ /  __/ /__/ /_/ /_/ / /    
                                     /_/  \__,_/____/\__/\__,_/     /_/   \___/\__, /_/\____/_/ /_/___  /_/_/ /_/____/ .___/\___/\___/\__/\____/_/     
                                                                             _/____/ __ \ |__ \  / __ \ / __ \      /_/                                
                                                                            | | / / / / / __/ / / / / // / / /                                         
                                                                            | |/ / /_/ / / __/_/ /_/ // /_/ /                                          
                                                                            |___/\____(_)____(_)____(_)____/                                           
                                                                                                                                                                                          
                                                                          Copyright (c) Matthew C. Mosior 2023       
  
[2023-12-09 23:09:17.706164190 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Querying BioMart for regions data.
[2023-12-09 23:09:17.706802837 EST] || LogInfo || ThreadId 4   || runQueryBioMart                                  || Generating BioMart compatible XML.
[2023-12-09 23:09:17.714542223 EST] || LogInfo || ThreadId 4   || runQueryBioMart                                  || Querying and downloading region data from BioMart via HTTP request.
[2023-12-09 23:09:19.972331939 EST] || LogInfo || ThreadId 4   || runQueryBioMart                                  || Successfully queried and returned region data from BioMart via HTTP request.
[2023-12-09 23:09:19.972426483 EST] || LogInfo || ThreadId 4   || runQueryBioMart                                  || Writing BioMart region data to file biomartresult.txt in output directory.
[2023-12-09 23:09:19.980867625 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Determining whether each variant is within its respective gene's TSS.
[2023-12-09 23:09:19.980920953 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Massaging TSS determination data into more usable format.
[2023-12-09 23:09:19.980954905 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Calculating the reverse complement of each user defined ambiguity code.
[2023-12-09 23:09:19.980981424 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Creating list of tuples to define directionality of each forward strand ambiguity code.
[2023-12-09 23:09:19.981005238 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Creating list of tuples to define directionality of each reverse strand ambiguity code.
[2023-12-09 23:09:19.981029022 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Generating all possible ambiguity code strings using SMT solver.
[2023-12-09 23:09:20.055093394 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Preparing ambiguity code strings to determine whether each lies within its respective TSS.
[2023-12-09 23:09:20.055156070 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Creating initial striped resource pool.
[2023-12-09 23:09:20.055233853 EST] || LogInfo || ThreadId 19  || fastaRegionInspect                               || Linearly processing all regions data.
[2023-12-09 23:09:20.055565934 EST] || LogInfo || ThreadId 20  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene BTG1.
[2023-12-09 23:09:20.055641995 EST] || LogInfo || ThreadId 20  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:20.055613662 EST] || LogInfo || ThreadId 21  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene MYH10.
[2023-12-09 23:09:20.055693139 EST] || LogInfo || ThreadId 21  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:20.055717985 EST] || LogInfo || ThreadId 22  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene ARID5B.
[2023-12-09 23:09:20.055743151 EST] || LogInfo || ThreadId 22  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for ARID5B.
[2023-12-09 23:09:20.056239896 EST] || LogInfo || ThreadId 22  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for ARID5B.
[2023-12-09 23:09:20.055775821 EST] || LogInfo || ThreadId 23  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene EBF1.
[2023-12-09 23:09:20.056303673 EST] || LogInfo || ThreadId 24  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene H1-5.
[2023-12-09 23:09:20.059131613 EST] || LogInfo || ThreadId 25  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene H1-2.
[2023-12-09 23:09:20.079130276 EST] || LogInfo || ThreadId 27  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene SGK1.
[2023-12-09 23:09:20.099131714 EST] || LogInfo || ThreadId 28  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene IL4R.
[2023-12-09 23:09:20.099147953 EST] || LogInfo || ThreadId 29  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene AFF3.
[2023-12-09 23:09:20.119436991 EST] || LogInfo || ThreadId 30  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene CXCR4.
[2023-12-09 23:09:20.119455916 EST] || LogInfo || ThreadId 23  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:20.139138957 EST] || LogInfo || ThreadId 31  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene GNA13.
[2023-12-09 23:09:20.159131328 EST] || LogInfo || ThreadId 32  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene IGLL5.
[2023-12-09 23:09:20.179138906 EST] || LogInfo || ThreadId 33  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene BCL7A.
[2023-12-09 23:09:20.199140313 EST] || LogInfo || ThreadId 34  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene CD83.
[2023-12-09 23:09:20.219139056 EST] || LogInfo || ThreadId 35  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene H2AC16.
[2023-12-09 23:09:20.239138650 EST] || LogInfo || ThreadId 36  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene TP53.
[2023-12-09 23:09:20.259140428 EST] || LogInfo || ThreadId 37  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene SOCS1.
[2023-12-09 23:09:20.339394025 EST] || LogInfo || ThreadId 24  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:20.459140217 EST] || LogInfo || ThreadId 25  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:20.579132060 EST] || LogInfo || ThreadId 27  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:20.659135816 EST] || LogInfo || ThreadId 28  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for IL4R.
[2023-12-09 23:09:20.699135034 EST] || LogInfo || ThreadId 29  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:20.779130625 EST] || LogInfo || ThreadId 30  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:20.839139606 EST] || LogInfo || ThreadId 31  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:20.879140288 EST] || LogInfo || ThreadId 32  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for IGLL5.
[2023-12-09 23:09:20.919139937 EST] || LogInfo || ThreadId 33  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for BCL7A.
[2023-12-09 23:09:20.959140677 EST] || LogInfo || ThreadId 34  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for CD83.
[2023-12-09 23:09:20.999140216 EST] || LogInfo || ThreadId 35  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for H2AC16.
[2023-12-09 23:09:21.039140529 EST] || LogInfo || ThreadId 36  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:21.079138827 EST] || LogInfo || ThreadId 37  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code WRC: WRC strand orientation is -1.
[2023-12-09 23:09:21.306975892 EST] || LogInfo || ThreadId 22  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for ARID5B.
[2023-12-09 23:09:21.307049638 EST] || LogInfo || ThreadId 28  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for IL4R.
[2023-12-09 23:09:21.307109989 EST] || LogInfo || ThreadId 32  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for IGLL5.
[2023-12-09 23:09:21.307140886 EST] || LogInfo || ThreadId 33  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for BCL7A.
[2023-12-09 23:09:21.307179658 EST] || LogInfo || ThreadId 34  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for CD83.
[2023-12-09 23:09:21.307215584 EST] || LogInfo || ThreadId 35  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for H2AC16.
[2023-12-09 23:09:21.307244207 EST] || LogInfo || ThreadId 22  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AGC.
[2023-12-09 23:09:21.659156727 EST] || LogInfo || ThreadId 22  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TGC.
[2023-12-09 23:09:21.659189598 EST] || LogInfo || ThreadId 22  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AAC.
[2023-12-09 23:09:21.659213291 EST] || LogInfo || ThreadId 22  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TAC.
[2023-12-09 23:09:22.883293704 EST] || LogInfo || ThreadId 34  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for CD83.
[2023-12-09 23:09:22.883336614 EST] || LogInfo || ThreadId 34  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AGC.
[2023-12-09 23:09:22.883353495 EST] || LogInfo || ThreadId 34  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TGC.
[2023-12-09 23:09:22.883369224 EST] || LogInfo || ThreadId 34  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AAC.
[2023-12-09 23:09:22.883383841 EST] || LogInfo || ThreadId 34  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TAC.
[2023-12-09 23:09:23.398343318 EST] || LogInfo || ThreadId 32  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for IGLL5.
[2023-12-09 23:09:23.398375357 EST] || LogInfo || ThreadId 32  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AGC.
[2023-12-09 23:09:23.398391396 EST] || LogInfo || ThreadId 32  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TGC.
[2023-12-09 23:09:23.398406655 EST] || LogInfo || ThreadId 32  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AAC.
[2023-12-09 23:09:23.398420490 EST] || LogInfo || ThreadId 32  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TAC.
[2023-12-09 23:09:23.515199757 EST] || LogInfo || ThreadId 28  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for IL4R.
[2023-12-09 23:09:23.515230654 EST] || LogInfo || ThreadId 28  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AGC.
[2023-12-09 23:09:23.515246313 EST] || LogInfo || ThreadId 28  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TGC.
[2023-12-09 23:09:23.516736408 EST] || LogInfo || ThreadId 28  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AAC.
[2023-12-09 23:09:23.516762286 EST] || LogInfo || ThreadId 28  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TAC.
[2023-12-09 23:09:23.754226286 EST] || LogInfo || ThreadId 35  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for H2AC16.
[2023-12-09 23:09:23.754255019 EST] || LogInfo || ThreadId 35  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AGC.
[2023-12-09 23:09:23.754283382 EST] || LogInfo || ThreadId 35  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TGC.
[2023-12-09 23:09:23.754301565 EST] || LogInfo || ThreadId 35  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AAC.
[2023-12-09 23:09:23.754325589 EST] || LogInfo || ThreadId 35  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TAC.
[2023-12-09 23:09:25.697450342 EST] || LogInfo || ThreadId 33  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for BCL7A.
[2023-12-09 23:09:25.697528466 EST] || LogInfo || ThreadId 33  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AGC.
[2023-12-09 23:09:25.697548964 EST] || LogInfo || ThreadId 33  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TGC.
[2023-12-09 23:09:25.697564773 EST] || LogInfo || ThreadId 33  || subStrLocationsSmallForward                      || Processing mapped ambiguity code AAC.
[2023-12-09 23:09:25.697579480 EST] || LogInfo || ThreadId 33  || subStrLocationsSmallForward                      || Processing mapped ambiguity code TAC.
[2023-12-09 23:09:25.697624182 EST] || LogInfo || ThreadId 39  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene BTG1.
[2023-12-09 23:09:25.697645081 EST] || LogInfo || ThreadId 39  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for BTG1.
[2023-12-09 23:09:25.697692668 EST] || LogInfo || ThreadId 39  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for BTG1.
[2023-12-09 23:09:25.699117664 EST] || LogInfo || ThreadId 40  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene MYH10.
[2023-12-09 23:09:25.699143191 EST] || LogInfo || ThreadId 40  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for MYH10.
[2023-12-09 23:09:25.699181932 EST] || LogInfo || ThreadId 40  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for MYH10.
[2023-12-09 23:09:25.719114479 EST] || LogInfo || ThreadId 41  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene ARID5B.
[2023-12-09 23:09:25.719143473 EST] || LogInfo || ThreadId 41  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code GYW: GYW strand orientation is 1.
[2023-12-09 23:09:25.719161466 EST] || LogInfo || ThreadId 42  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene EBF1.
[2023-12-09 23:09:25.719178688 EST] || LogInfo || ThreadId 42  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for EBF1.
[2023-12-09 23:09:25.719271629 EST] || LogInfo || ThreadId 42  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for EBF1.
[2023-12-09 23:09:25.719202021 EST] || LogInfo || ThreadId 43  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene H1-5.
[2023-12-09 23:09:25.719932296 EST] || LogInfo || ThreadId 43  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for H1-5.
[2023-12-09 23:09:25.721589810 EST] || LogInfo || ThreadId 43  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for H1-5.
[2023-12-09 23:09:25.839114836 EST] || LogInfo || ThreadId 49  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene H1-2.
[2023-12-09 23:09:25.839170578 EST] || LogInfo || ThreadId 49  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for H1-2.
[2023-12-09 23:09:25.839267767 EST] || LogInfo || ThreadId 49  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for H1-2.
[2023-12-09 23:09:25.839202808 EST] || LogInfo || ThreadId 50  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene SGK1.
[2023-12-09 23:09:25.839326746 EST] || LogInfo || ThreadId 50  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for SGK1.
[2023-12-09 23:09:25.839826737 EST] || LogInfo || ThreadId 50  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for SGK1.
[2023-12-09 23:09:25.999118761 EST] || LogInfo || ThreadId 53  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene IL4R.
[2023-12-09 23:09:25.999147634 EST] || LogInfo || ThreadId 53  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code GYW: GYW strand orientation is 1.
[2023-12-09 23:09:25.999167961 EST] || LogInfo || ThreadId 54  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene AFF3.
[2023-12-09 23:09:25.999186115 EST] || LogInfo || ThreadId 54  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for AFF3.
[2023-12-09 23:09:25.999247708 EST] || LogInfo || ThreadId 54  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for AFF3.
[2023-12-09 23:09:26.139117143 EST] || LogInfo || ThreadId 55  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene CXCR4.
[2023-12-09 23:09:26.139145876 EST] || LogInfo || ThreadId 55  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for CXCR4.
[2023-12-09 23:09:26.139188945 EST] || LogInfo || ThreadId 55  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for CXCR4.
[2023-12-09 23:09:26.159115983 EST] || LogInfo || ThreadId 56  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene GNA13.
[2023-12-09 23:09:26.159152250 EST] || LogInfo || ThreadId 56  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for GNA13.
[2023-12-09 23:09:26.159261461 EST] || LogInfo || ThreadId 56  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for GNA13.
[2023-12-09 23:09:26.159175272 EST] || LogInfo || ThreadId 57  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene IGLL5.
[2023-12-09 23:09:26.160056856 EST] || LogInfo || ThreadId 57  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code GYW: GYW strand orientation is 1.
[2023-12-09 23:09:26.359115254 EST] || LogInfo || ThreadId 62  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene BCL7A.
[2023-12-09 23:09:26.359173221 EST] || LogInfo || ThreadId 62  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code GYW: GYW strand orientation is 1.
[2023-12-09 23:09:26.359198558 EST] || LogInfo || ThreadId 63  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene CD83.
[2023-12-09 23:09:26.359220027 EST] || LogInfo || ThreadId 63  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code GYW: GYW strand orientation is 1.
[2023-12-09 23:09:26.359240786 EST] || LogInfo || ThreadId 64  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene H2AC16.
[2023-12-09 23:09:26.359260592 EST] || LogInfo || ThreadId 64  || ambiguityCodesWithinRegionCheckSmall             || Could not process region data associated with current ambiguity code GYW: GYW strand orientation is 1.
[2023-12-09 23:09:26.359280018 EST] || LogInfo || ThreadId 65  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene TP53.
[2023-12-09 23:09:26.359299674 EST] || LogInfo || ThreadId 65  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for TP53.
[2023-12-09 23:09:26.539126485 EST] || LogInfo || ThreadId 66  || ambiguityCodesWithinRegionCheckSmall             || Processing region data associated with gene SOCS1.
[2023-12-09 23:09:26.539186766 EST] || LogInfo || ThreadId 66  || ambiguityCodesWithinRegionCheckSmall             || Reading fasta index file for SOCS1.
[2023-12-09 23:09:26.539261083 EST] || LogInfo || ThreadId 65  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for TP53.
[2023-12-09 23:09:26.539345368 EST] || LogInfo || ThreadId 66  || ambiguityCodesWithinRegionCheckSmall             || Reading in fasta file for SOCS1.
[2023-12-09 23:09:26.973299648 EST] || LogInfo || ThreadId 40  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for MYH10.
[2023-12-09 23:09:27.639154329 EST] || LogInfo || ThreadId 40  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:27.639213027 EST] || LogInfo || ThreadId 40  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:27.639237913 EST] || LogInfo || ThreadId 40  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:27.639261085 EST] || LogInfo || ThreadId 40  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:28.271554270 EST] || LogInfo || ThreadId 65  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for TP53.
[2023-12-09 23:09:28.271602509 EST] || LogInfo || ThreadId 65  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:28.271621534 EST] || LogInfo || ThreadId 65  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:28.271638154 EST] || LogInfo || ThreadId 65  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:28.271654184 EST] || LogInfo || ThreadId 65  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:29.510626754 EST] || LogInfo || ThreadId 66  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for SOCS1.
[2023-12-09 23:09:29.510683137 EST] || LogInfo || ThreadId 66  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:29.511787352 EST] || LogInfo || ThreadId 66  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:29.511810074 EST] || LogInfo || ThreadId 66  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:29.511829099 EST] || LogInfo || ThreadId 66  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:30.409990602 EST] || LogInfo || ThreadId 49  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for H1-2.
[2023-12-09 23:09:30.410051364 EST] || LogInfo || ThreadId 49  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:30.410069087 EST] || LogInfo || ThreadId 49  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:30.410084435 EST] || LogInfo || ThreadId 49  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:30.410098741 EST] || LogInfo || ThreadId 49  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:30.633488908 EST] || LogInfo || ThreadId 43  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for H1-5.
[2023-12-09 23:09:30.633521999 EST] || LogInfo || ThreadId 43  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:30.633544240 EST] || LogInfo || ThreadId 43  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:30.633560921 EST] || LogInfo || ThreadId 43  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:30.633577101 EST] || LogInfo || ThreadId 43  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:35.614113264 EST] || LogInfo || ThreadId 56  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for GNA13.
[2023-12-09 23:09:35.614171040 EST] || LogInfo || ThreadId 56  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:35.614643681 EST] || LogInfo || ThreadId 56  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:35.614665502 EST] || LogInfo || ThreadId 56  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:35.614682653 EST] || LogInfo || ThreadId 56  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:38.132036432 EST] || LogInfo || ThreadId 39  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for BTG1.
[2023-12-09 23:09:38.132089360 EST] || LogInfo || ThreadId 39  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:38.133120951 EST] || LogInfo || ThreadId 39  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:38.133145967 EST] || LogInfo || ThreadId 39  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:38.133166535 EST] || LogInfo || ThreadId 39  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:39.306229117 EST] || LogInfo || ThreadId 54  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for AFF3.
[2023-12-09 23:09:39.306284479 EST] || LogInfo || ThreadId 54  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:39.306303474 EST] || LogInfo || ThreadId 54  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:39.306317400 EST] || LogInfo || ThreadId 54  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:39.306330584 EST] || LogInfo || ThreadId 54  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:41.308159091 EST] || LogInfo || ThreadId 50  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for SGK1.
[2023-12-09 23:09:41.308217258 EST] || LogInfo || ThreadId 50  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:41.308235782 EST] || LogInfo || ThreadId 50  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:41.308251431 EST] || LogInfo || ThreadId 50  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:41.308266048 EST] || LogInfo || ThreadId 50  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:41.506848358 EST] || LogInfo || ThreadId 55  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for CXCR4.
[2023-12-09 23:09:41.506875788 EST] || LogInfo || ThreadId 55  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:41.506891938 EST] || LogInfo || ThreadId 55  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:41.506907998 EST] || LogInfo || ThreadId 55  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:41.506923326 EST] || LogInfo || ThreadId 55  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:41.955991211 EST] || LogInfo || ThreadId 42  || ambiguityCodesWithinRegionCheckSmall             || Grabbing all mapped ambiguity string locations for EBF1.
[2023-12-09 23:09:41.956031806 EST] || LogInfo || ThreadId 42  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCT.
[2023-12-09 23:09:41.956051702 EST] || LogInfo || ThreadId 42  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTA.
[2023-12-09 23:09:41.956066670 EST] || LogInfo || ThreadId 42  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GCA.
[2023-12-09 23:09:41.956081387 EST] || LogInfo || ThreadId 42  || subStrLocationsSmallReverse                      || Processing mapped ambiguity code GTT.
[2023-12-09 23:09:41.956355753 EST] || LogInfo || ThreadId 69  || fastaRegionInspect                               || Processing all variant data.
[2023-12-09 23:09:41.956536396 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Preparing variants for final analysis.
[2023-12-09 23:09:41.981279050 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Preparing to produce output CSV files.
[2023-12-09 23:09:41.981393711 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Producing output CSV files.
[2023-12-09 23:09:41.986557731 EST] || LogInfo || ThreadId 4   || fastaRegionInspect                               || Shutting down fasta-region-inspector v0.2.0.0.
```

## Dependency Visualization via Stack Dot
The following is a visualization of the dependency graph for **FRI**:


![](dependencies.png)
