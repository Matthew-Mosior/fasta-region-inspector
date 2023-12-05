{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE KindSignatures    #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE RankNTypes        #-}
{-# LANGUAGE Strict            #-}
{-# LANGUAGE TypeOperators     #-}

module Inspect where

import Amalgamate
import AmbiguityCodes
import CmdOpts
import Logging
import QueryBioMart
import Region
import Utility
import Variant
import Types

import Linear.RegionsLinear
import Linear.VariantLinear
import Linear.UtilityLinear

import           Data.Aeson.Types
import           Data.ByteString as DB
import           Data.ByteString.Char8 as DBC
import           Data.ByteString.Lazy as DBL
import           Data.ByteString.Search.DFA as DBSDFA
import           Data.Char as DC
import           Data.Functor as DF
import           Data.List as DL
import           Data.List.Split as DLS
import           Data.Ord as DO
import           Data.SBV as DSBV
import qualified Data.SBV.String as DSBVS
import           Data.SBV.RegExp as DSBVRE
import qualified Data.Text as DText
import           Data.Time as DTime
import           Data.Traversable as DT
import           Effectful
import           Effectful.Ki
import           GHC.Stack                            (callStack,getCallStack,HasCallStack)
import           UnliftIO.Pool                        (mkDefaultPoolConfig,newPool,withResource)

fastaRegionInspect :: forall {es :: [Effect]} {b}.
                      ( StructuredConcurrency :> es
                      , IOE :> es
                      , HasCallStack
                      )
                   => FRIConfig
                   -> Eff es ()
fastaRegionInspect config = do
  --Get calling function name.
  let ((callingfunction,_):_) = getCallStack callStack
  --Query biomart for Region data.
  _ <- liftIO $ showPrettyLog LogInfo
                              (maxnumberconcthreads config)
                              callingfunction
                              "Querying BioMart for regions data."
  biomartregiondata <- runQueryBioMart config
  --Determine whether each variant is within
  --the TSS of the genes each variant is located in.
  _ <- liftIO $ showPrettyLog LogInfo
                              (maxnumberconcthreads config)
                              callingfunction
                              "Determining whether each variant is within its respective gene's TSS."
  let withintss = variantWithinRegionCheck config
                                           biomartregiondata
  --Prepare withintss, ambiguitycodeswithintss, and variantsinambiguitycodesandtss for printing.
  _ <- liftIO $ showPrettyLog LogInfo
                              (maxnumberconcthreads config)
                              callingfunction
                              "Massaging TSS determination data into more usable format."
  let printreadywithintss = prepareWithinTSS withintss
  --Grab the reverse complement of the 
  --user defined ambiguity codes.
  _ <- liftIO $ showPrettyLog LogInfo
                              (maxnumberconcthreads config)
                              callingfunction
                              "Calculating the reverse complement of each user defined ambiguity code."
  let ambiguitycodesreversecomplements = ambiguityCodesReverseComplement config
  --Create list of tuples defining directionality of ambiguitycodes.
  _ <- liftIO $ showPrettyLog LogInfo
                              (maxnumberconcthreads config)
                              callingfunction
                              "Creating list of tuples to define directionality of each forward strand ambiguity code."
  let ambiguitycodesfinaltuple = DL.map (\x -> (x,"1"))
                                 (DL.map (DText.unpack) (ambiguitycodes config))
  --Create list of tuples defining directionality of ambiguitycodesreversecomplements.
  _ <- liftIO $ showPrettyLog LogInfo
                              (maxnumberconcthreads config)
                              callingfunction
                              "Creating list of tuples to define directionality of each reverse strand ambiguity code."
  let ambiguitycodesreversecomplementstuple = DL.map (\x -> (x,"-1"))
                                              ambiguitycodesreversecomplements 
  --Grab all possible strings created from each ambiguity codes.
  _ <- liftIO $ showPrettyLog LogInfo
                              (maxnumberconcthreads config)
                              callingfunction
                              "Generating all possible ambiguity code strings using SMT solver."
  allmappedambiguitystrs <- liftIO $ allStrGeneration ( DL.map DText.unpack (ambiguitycodes config) ++
                                                        ambiguitycodesreversecomplements
                                                      ) 
  --Prepare allmappedambiguitystrs for ambiguityCodesWithinRegionCheck.
  _ <- liftIO $ showPrettyLog LogInfo
                              (maxnumberconcthreads config)
                              callingfunction
                              "Preparing ambiguity code strings to determine whether each lies within its respective TSS."
  let allmappedambiguitystrstuple = stringToTuple allmappedambiguitystrs
  --Initialize a resource pool (for controlling threading).
  _ <- liftIO $ showPrettyLog LogInfo
                              (maxnumberconcthreads config)
                              callingfunction
                              "Creating initial striped resource pool." 
  resourcepoolconfig <- liftIO $ mkDefaultPoolConfig (pure ())
                                                     (\_ -> pure ())
                                                     (0.5 :: Double)
                                                     (maxnumberconcthreads config)
  newresourcepool    <- liftIO $ newPool resourcepoolconfig
  --Fork a scoped thread via effecful-ki
  --for variant and region analysis.
  withResource newresourcepool $ \_ -> do
    scoped $ \scope -> do
      (ambiguitycodeswithintss,analysisreadyambiguitycodeswithintss) <- do
        regions <- fork scope (do --Utilize linear resources to open the input
                                  --fasta file and read in the sequences associated with the
                                  --regions of interest.
                                  _ <- liftIO $ showPrettyLog LogInfo
                                                              (maxnumberconcthreads config)
                                                              callingfunction
                                                              "Linearly processing all regions data."
                                  regionsLinear (tsswindowsize config)
                                                config
                                                (ambiguitycodesfinaltuple ++ ambiguitycodesreversecomplementstuple)
                                                allmappedambiguitystrstuple
                                                biomartregiondata
                              )
        atomically $ await regions 
      variantsf <- do variants <- fork scope (do --Utilize linear resources to open the input
                                                 --fasta file and read in the sequence associated with
                                                 --the variants sequence description.
                                                 _ <- liftIO $ showPrettyLog LogInfo
                                                                             (maxnumberconcthreads config)
                                                                             callingfunction
                                                                             "Processing all variant data."
                                                 variantLinear withintss
                                                               analysisreadyambiguitycodeswithintss
                                             )
                      --Wait for vaiants thread to terminate.
                      atomically $ await variants
      --Prepare withintss, ambiguitycodeswithintss, and variantsinambiguitycodesandtss for printing.
      _ <- liftIO $ showPrettyLog LogInfo
                                  (maxnumberconcthreads config)
                                  callingfunction
                                  "Preparing variants for final analysis."
      let printreadywithintss = prepareWithinTSS withintss
      let printreadyambiguitycodeswithintss = DL.concat
                                              (DL.map
                                              (DL.map
                                              (\xs -> (DL.take 6 xs)
                                              ++ [DL.intercalate "," (DL.drop 6 xs)]))
                                              (prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss))
      --Prepare final amalgamated file.
      let finalvariantfile = amalgamateFinalVariantData printreadywithintss
                                                        variantsf
      --Prepare final print ready files with headers.
      _ <- liftIO $ showPrettyLog LogInfo
                                  (maxnumberconcthreads config)
                                  callingfunction
                                  "Preparing to produce output CSV files."
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
      _ <- liftIO $ showPrettyLog LogInfo
                                  (maxnumberconcthreads config)
                                  callingfunction
                                  "Producing output CSV files."
      if | writeambiguitycodes config
         -> do liftIO $ writeCSV config
                                 ambiguitycodeswithintsscsv
                                 "ambiguity_codes.csv"
               liftIO $ writeCSV config
                                 finalvariantfilecsv
                                 "variants_in_ambiguity_codes.csv"
               --Shut down FRI.
               liftIO $ showPrettyLog LogInfo
                                      (maxnumberconcthreads config)
                                      callingfunction
                                      "Shutting down fasta-region-inspector v0.2.0.0."
         | otherwise
         -> do liftIO $ writeCSV config
                                 finalvariantfilecsv 
                                 "variants_in_ambiguity_codes.csv"
               --Shut down FRI.
               liftIO $ showPrettyLog LogInfo
                                      (maxnumberconcthreads config)
                                      callingfunction
                                      "Shutting down fasta-region-inspector v0.2.0.0."
