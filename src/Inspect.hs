{-=Fasta-Region-Inspector (FRI): A Somatic=-}
{-=Hypermutation Analysis Tool.=-}
{-=Author: Matthew Mosior=-}
{-=Synposis: This Haskell script will=-}
{-=process command line arguments to FRI.=-}


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

import Control.Parallel.Strategies as CPS
import Data.Attoparsec.ByteString as DAB
import Data.ByteString as DB
import Data.ByteString.Char8 as DBC
import Data.ByteString.Lazy as DBL
import Data.ByteString.Search.DFA as DBSDFA
import Data.Compact as DCompact
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
import ELynx.Sequence.Import.Fasta as EISF
import qualified ELynx.Alphabet.Alphabet as EDAA
import qualified ELynx.Alphabet.Character as EDAC
import qualified ELynx.Character.Character as EDCC
import qualified ELynx.Character.NucleotideI as EDCN
import qualified ELynx.Sequence.Sequence as EDSS
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
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Starting up Fasta Region Inspector v0.1.0.0 ...")
  --Query BioMart for Region data.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Query BioMart for regions data ...")
  biomartregiondata <- runQueryBioMart config
  --Determine whether each variant is within
  --the TSS of the genes each variant is located in.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Determining whether each variant is within the its respective genes TSS ...")
  let withintss = variantWithinRegionCheck config
                                           biomartregiondata
  --Grab the reverse complement of the 
  --user defined ambiguity codes.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Calculating the reverse complement of each user defined ambiguity code ...")
  let ambiguitycodesreversecomplements = ambiguityCodesReverseComplement config
  --Create list of tuples defining directionality of ambiguitycodes.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Creating list of tuples to define directionality of each forward strand ambiguity code ...")
  let ambiguitycodesfinaltuple = DL.map (\x -> (x,"1"))
                                 (DL.map (DText.unpack) (ambiguitycodes config))
  --Create list of tuples defining directionality of ambiguitycodesreversecomplements.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Creating list of tuples to define directionality of each reverse strand ambiguity code ...")
  let ambiguitycodesreversecomplementstuple = DL.map (\x -> (x,"-1"))
                                              ambiguitycodesreversecomplements 
  --Grab all possible strings created from each ambiguity codes.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Generating all possible ambiguity codes strings using SMT solver ...")
  allmappedambiguitystrs <- allStrGeneration ( DL.map (DText.unpack) (ambiguitycodes config) ++
                                               ambiguitycodesreversecomplements
                                             ) 
  --Prepare allmappedambiguitystrs for ambiguityCodesWithinRegionCheck.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Preparing ambiguity code strings to determine whether each lies within its respective TSS ...")
  let allmappedambiguitystrstuple = stringToTuple allmappedambiguitystrs
  --Read fasta file into a strict ByteString.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Reading fasta file into strict ByteString ...")
  inputfastafile <- DB.readFile  $
                    DText.unpack $
                    YamlParser.fasta config
  --Parse input fasta file to [Sequence].
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Parsing fasta file into [Sequence] ...")
  let pfastafile = parseOnly (EISF.fasta EDAA.DNAI <* endOfInput)
                             inputfastafile
  case pfastafile of
    Left  err    -> do currenttandd <- DTime.getZonedTime
                       _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                              ++ "Could not parse fasta file into [Sequence]: "
                                              ++ err)
                       currenttandd <- DTime.getZonedTime
                       _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                              ++ "Shutting down Fasta Region Inspector v0.1.0.0 ...")
                       SX.exitWith (SX.ExitFailure 1)
    Right cfasta -> do --Extract every sequences name and characters from cfasta.
                       currenttandd <- DTime.getZonedTime
                       _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                              ++ "Extract every sequences name and characters from cfasta ...")
                       let namesandcharacters = extractNameAndCharacters cfasta
                       --Put namesandcharacters into compact region.
                       currenttandd <- DTime.getZonedTime
                       _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                              ++ "Putting namesandcharacters into compact region ...")
                       namesandcharactersc <- DCompact.compact namesandcharacters
                       --Determine whether there are ambiguity codes strings
                       --present within the TSS of each region.
                       currenttandd <- DTime.getZonedTime
                       _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                              ++ "Determing whether each ambiguity code string lies within its respective TSS ...") 
                       if | ignorestrandedness config
                          -> do ambiguitycodeswithintss <- ambiguityCodesWithinRegionCheckIgnoreStrand
                                                           (tsswindowsize config)
                                                           (DCompact.getCompact namesandcharactersc)
                                                           (ambiguitycodesfinaltuple ++ ambiguitycodesreversecomplementstuple) 
                                                           allmappedambiguitystrstuple
                                                           biomartregiondata  
                                --Prepare ambiguitycodeswithintss for printing.
                                currenttandd <- DTime.getZonedTime
                                _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                       ++ "Preparing ambiguity code strings final analysis ...")
                                let analysisreadyambiguitycodeswithintss = prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss 
                                --Determine whether there are variants present
                                --within ambiguity codes within corresponding regions.
                                let printreadyvariantsinambiguitycodesandtss = variantsWithinAmbiguityCodesAndTSS withintss
                                                                                                                  analysisreadyambiguitycodeswithintss
                                --Prepare withintss, ambiguitycodeswithintss, and variantsinambiguitycodesandtss for printing.
                                currenttandd <- DTime.getZonedTime
                                _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                       ++ "Preparing variants final analysis ...")
                                let printreadywithintss = prepareWithinTSS withintss
                                let printreadyambiguitycodeswithintss = DL.concat 
                                                                        (DL.map 
                                                                        (DL.map 
                                                                        (\xs -> (DL.take 6 xs) 
                                                                        ++ [DL.intercalate "," (DL.drop 6 xs)])) 
                                                                        (prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss))
                                --Prepare final amalgamated file.
                                let finalvariantfile = amalgamateFinalVariantData printreadywithintss
                                                                                  printreadyvariantsinambiguitycodesandtss 
                                --Prepare final print ready files with headers.
                                currenttandd <- DTime.getZonedTime
                                _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                       ++ "Preparing to print output CSV files ...") 
                                let finalprintreadyambiguitycodeswithintss = [ "Ambiguity_Code"
                                                                             , "Mapped_Nucleotide_String"
                                                                             , "Chromosome"
                                                                             , "TSS"
                                                                             , "Strand"
                                                                             , "SYMBOL"
                                                                             , "Ambiguity_Code_String_Locations_Within_TSS"
                                                                             ] 
                                                                             : printreadyambiguitycodeswithintss
                                let ambiguitycodeswithintsscsv = toCSV finalprintreadyambiguitycodeswithintss
                                let printreadyfinalvariantfile = [ "Variant"
                                                                 , "Region"
                                                                 , "Variant_Within_Region"
                                                                 , "Ambiguity_Code"
                                                                 , "Mapped_Nucleotide_String"
                                                                 , "Ambiguity_Code_String_Locations_Within_TSS"
                                                                 ]
                                                                 : finalvariantfile  
                                let finalvariantfilecsv = toCSV printreadyfinalvariantfile
                                --Print withintss, ambiguitycodeswithintss, and variantsinambiguitycodesandtss to files.
                                currenttandd <- DTime.getZonedTime
                                _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                       ++ "Printing output CSV files ...")
                                if | writeambiguitycodes config
                                   -> do writeCSV config
                                                  ambiguitycodeswithintsscsv
                                                  "ambiguity_codes.csv"
                                         writeCSV config
                                                  finalvariantfilecsv 
                                                  "variants_in_ambiguity_codes.csv"
                                         --Shut down FRI.
                                         currenttandd <- DTime.getZonedTime
                                         SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                           ++ "Shutting down Fasta Region Inspector v0.1.0.0 ...")
                                   | otherwise
                                   -> do writeCSV config
                                                  finalvariantfilecsv 
                                                  "variants_in_ambiguity_codes.csv"
                                         --Shut down FRI.
                                         currenttandd <- DTime.getZonedTime
                                         SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                           ++ "Shutting down Fasta Region Inspector v0.1.0.0 ...")
                          | otherwise
                          -> do ambiguitycodeswithintss <- ambiguityCodesWithinRegionCheck
                                                           (tsswindowsize config)
                                                           (DCompact.getCompact namesandcharactersc)
                                                           (ambiguitycodesfinaltuple ++ ambiguitycodesreversecomplementstuple)
                                                           allmappedambiguitystrstuple
                                                           biomartregiondata
                                --Prepare ambiguitycodeswithintss for printing.
                                currenttandd <- DTime.getZonedTime
                                _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                       ++ "Preparing ambiguity code strings final analysis ...")
                                let analysisreadyambiguitycodeswithintss = prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss
                                --Determine whether there are variants present
                                --within ambiguity codes within corresponding regions.
                                let printreadyvariantsinambiguitycodesandtss = variantsWithinAmbiguityCodesAndTSS withintss
                                                                                                                  analysisreadyambiguitycodeswithintss
                                --Prepare withintss, ambiguitycodeswithintss, and variantsinambiguitycodesandtss for printing.
                                currenttandd <- DTime.getZonedTime
                                _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                       ++ "Preparing variants final analysis ...")
                                let printreadywithintss = prepareWithinTSS withintss
                                let printreadyambiguitycodeswithintss = DL.concat
                                                                        (DL.map
                                                                        (DL.map
                                                                        (\xs -> (DL.take 6 xs)
                                                                        ++ [DL.intercalate "," (DL.drop 6 xs)]))
                                                                        (prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss))
                                --Prepare final amalgamated file.
                                let finalvariantfile = amalgamateFinalVariantData printreadywithintss
                                                                                  printreadyvariantsinambiguitycodesandtss 
                                --Prepare final print ready files with headers.
                                currenttandd <- DTime.getZonedTime
                                _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                       ++ "Preparing to print output CSV files ...")
                                let finalprintreadyambiguitycodeswithintss = [ "Ambiguity_Code"
                                                                             , "Mapped_Nucleotide_String"
                                                                             , "Chromosome"
                                                                             , "TSS"
                                                                             , "Strand"
                                                                             , "SYMBOL"
                                                                             , "Ambiguity_Code_String_Locations_Within_TSS"
                                                                             ]
                                                                             : printreadyambiguitycodeswithintss
                                let ambiguitycodeswithintsscsv = toCSV finalprintreadyambiguitycodeswithintss
                                let printreadyfinalvariantfile = [ "Variant"
                                                                 , "Region"
                                                                 , "Variant_Within_Region"
                                                                 , "Ambiguity_Code"
                                                                 , "Mapped_Nucleotide_String"
                                                                 , "Ambiguity_Code_String_Locations_Within_TSS"
                                                                 ]
                                                                 : finalvariantfile 
                                --let variantsinambiguitycodesandtsscsv = toCSV finalprintreadyvariantsinambiguitycodesandtss
                                let finalvariantfilecsv = toCSV printreadyfinalvariantfile 
                                --Print withintss, ambiguitycodeswithintss, and variantsinambiguitycodesandtss to files.
                                currenttandd <- DTime.getZonedTime
                                _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                       ++ "Printing output CSV files ...")
                                if | writeambiguitycodes config
                                   -> do writeCSV config
                                                  ambiguitycodeswithintsscsv
                                                  "ambiguity_codes.csv"
                                         writeCSV config
                                                  finalvariantfilecsv
                                                  "variants_in_ambiguity_codes.csv"
                                         --Shut down FRI.
                                         currenttandd <- DTime.getZonedTime
                                         SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                           ++ "Shutting down Fasta Region Inspector v0.1.0.0 ...")
                                   | otherwise
                                   -> do writeCSV config
                                                  finalvariantfilecsv 
                                                  "variants_in_ambiguity_codes.csv"
                                         --Shut down FRI.
                                         currenttandd <- DTime.getZonedTime
                                         SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                           ++ "Shutting down Fasta Region Inspector v0.1.0.0 ...")

{--------------------}
