{-=Fasta-Region-Inspector (FRI): A Somatic=-}
{-=Hypermutation Analysis Tool.=-}
{-=Author: Matthew Mosior=-}
{-=Synposis: This Haskell script will=-}
{-=process command line arguments to FRI.=-}


{-Language extension.-}

{-# LANGUAGE Strict     #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE MultiWayIf #-}

{---------------------}


{-Module.-}

module Inspect where

{---------}


{-Import modules.-}

import CmdOpts
import Common
import QueryBioMart
import YamlParser

{-----------------}


{-Imports.-}

import Control.DeepSeq as CD
import Control.Parallel.Strategies as CPS
import Data.ByteString as DB
import Data.ByteString.Char8 as DBC
import Data.ByteString.Lazy as DBL
import Data.ByteString.Search.DFA as DBSDFA
import Data.Char as DC
import Data.Functor as DF
import Data.List as DL
import Data.List.Split as DLS
import Data.Ord as DO
import Data.SBV as DSBV
import qualified Data.SBV.String as DSBVS
import Data.SBV.RegExp as DSBVRE
import Data.Text as DText
import Data.Time as DTime
import Data.Traversable as DT
import System.Console.GetOpt as SCG
import System.Process as SP
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO
import System.IO.Temp as SIOT
import Text.Regex.TDFA as TRTDFA

{----------}


{-Inspect functions.-}

fastaRegionInspect :: FRIConfig -> IO ()
fastaRegionInspect config = do
  --Starting up Fasta Region Inspector v0.1.0.0.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Starting up Fasta Region Inspector v0.1.0.0 ...")
  --Query BioMart for Region data.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Query BioMart for regions data ...")
  biomartregiondata <- runQueryBioMart config
  --Determine whether each variant is within
  --the TSS of the genes each variant is located in.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Determining whether each variant is within the its respective genes TSS ...")
  let withintss = variantWithinRegionCheck config
                                           biomartregiondata
  --Grab the reverse complement of the 
  --user defined ambiguity codes.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Calculating the reverse complement of each user defined ambiguity code ...")
  let ambiguitycodesreversecomplements = ambiguityCodesReverseComplement config
  --Create list of tuples defining directionality of ambiguitycodes.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Creating list of tuples to define directionality of each forward strand ambiguity code ...")
  let ambiguitycodesfinaltuple = DL.map (\x -> (x,"1"))
                                 (DL.map (DText.unpack) (ambiguitycodes config))
  --Create list of tuples defining directionality of ambiguitycodesreversecomplements.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Creating list of tuples to define directionality of each reverse strand ambiguity code ...")
  let ambiguitycodesreversecomplementstuple = DL.map (\x -> (x,"-1"))
                                              ambiguitycodesreversecomplements 
  --Grab all possible strings created from each ambiguity codes.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Generating all possible ambiguity codes strings using SMT solver ...")
  allmappedambiguitystrs <- allStrGeneration ( DL.map (DText.unpack) (ambiguitycodes config) ++
                                               ambiguitycodesreversecomplements
                                             ) 
  --Prepare allmappedambiguitystrs for ambiguityCodesWithinRegionCheck.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Preparing ambiguity code strings to determine whether each lies within its respective TSS ...")
  let allmappedambiguitystrstuple = stringToTuple allmappedambiguitystrs
  --Determine whether there are ambiguity codes strings
  --present within the TSS of each region.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Determing whether each ambiguity code string lies within its respective TSS ...") 
  ambiguitycodeswithintss <- ambiguityCodesWithinRegionCheck
                             config 
                             (ambiguitycodesfinaltuple ++ ambiguitycodesreversecomplementstuple) 
                             allmappedambiguitystrstuple
                             biomartregiondata  
  --Prepare ambiguitycodeswithintss for printing.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Preparing ambiguity code strings final analysis ...")
  let analysisreadyambiguitycodeswithintss = prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss 
  --Determine whether there are variants present
  --within ambiguity codes within corresponding regions.
  let printreadyvariantsinambiguitycodesandtss = ambiguitycodeswithintss 
                                                 `CD.deepseq` 
                                                 variantsWithinAmbiguityCodesAndTSS withintss analysisreadyambiguitycodeswithintss
  --Prepare withintss, ambiguitycodeswithintss, and variantsinambiguitycodesandtss for printing.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Preparing variants final analysis ...")
  let printreadywithintss = prepareWithinTSS withintss
  let printreadyambiguitycodeswithintss = ambiguitycodeswithintss 
                                          `CD.deepseq` 
                                          DL.concat 
                                          (DL.map 
                                          (DL.map 
                                          (\xs -> (DL.take 6 xs) 
                                          ++ [DL.intercalate "," (DL.drop 6 xs)])) 
                                          (prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss)) 
  --Prepare final print ready files with headers.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Preparing to print output CSV files ...")
  let finalprintreadywithintss = printreadywithintss
                                 `CD.deepseq`
                                 [["Variant"
                                  ,"Region"
                                  ,"Variant_Within_Region"]] 
                                 ++ printreadywithintss
  let withintsscsv = finalprintreadywithintss
                     `CD.deepseq`
                     toCSV finalprintreadywithintss
  let finalprintreadyambiguitycodeswithintss =  printreadyambiguitycodeswithintss
                                                `CD.deepseq`
                                                [["Ambiguity_Code"
                                                ,"Mapped_Nucleotide_String"
                                                ,"Chromosome"
                                                ,"TSS"
                                                ,"Strand"
                                                ,"SYMBOL"
                                                ,"Ambiguity_Code_String_Locations_Within_TSS"]] 
                                               ++ printreadyambiguitycodeswithintss
  let ambiguitycodeswithintsscsv = finalprintreadyambiguitycodeswithintss
                                   `CD.deepseq`
                                   toCSV finalprintreadyambiguitycodeswithintss
  let finalprintreadyvariantsinambiguitycodesandtss = printreadyvariantsinambiguitycodesandtss
                                                      `CD.deepseq`
                                                      [["Variant"
                                                       ,"Region"
                                                       ,"Variant_Within_Region"
                                                       ,"Ambiguity_Code"
                                                       ,"Mapped_Nucleotide_String"
                                                       ,"Ambiguity_Code_String_Locations_Within_TSS"]] 
                                                      ++ printreadyvariantsinambiguitycodesandtss
  let variantsinambiguitycodesandtsscsv = finalprintreadyvariantsinambiguitycodesandtss
                                          `CD.deepseq`
                                          toCSV finalprintreadyvariantsinambiguitycodesandtss
  --Print withintss, ambiguitycodeswithintss, and variantsinambiguitycodesandtss to files.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Printing output CSV files ...") 
  finalprintreadywithintss `CD.deepseq` writeCSV config
                                                 withintsscsv
                                                 "variants.csv"
  finalprintreadyambiguitycodeswithintss `CD.deepseq` writeCSV config
                                                               ambiguitycodeswithintsscsv
                                                               "ambiguity_codes.csv"
  finalprintreadyvariantsinambiguitycodesandtss `CD.deepseq` writeCSV config
                                                                      variantsinambiguitycodesandtsscsv
                                                                      "variants_in_ambiguity_codes.csv"
  --Shut down FRI.
  currenttandd <- DTime.getZonedTime
  SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                    ++ "Shutting down Fasta Region Inspector v0.1.0.0 ...")

{--------------------}
