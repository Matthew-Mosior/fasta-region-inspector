{-=Fasta-Region-Inspector (FRI): A Somatic=-}
{-=Hypermutation Analysis Tool.=-}
{-=Author: Matthew Mosior=-}
{-=Synposis: This Haskell script describes=-}
{-=the configuration YAML file.=-}


{-Language extension.-}

{-# LANGUAGE Strict            #-}
{-# LANGUAGE StrictData        #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE RecordWildCards   #-}
{-# LANGUAGE DeriveGeneric     #-}

{---------------------}


{-Module.-}

module YamlParser where

{---------}


{-Imports.-}

import Data.Aeson as DAeson
import Data.HashMap.Lazy as DHL
import Data.List as DL
import Data.Maybe as DMaybe
import Data.Text as DText
import Data.Yaml as DYaml
import Control.Applicative as CA
import GHC.Generics 

{----------}


{-Custom YAML input file Datatype and related functions.-}

data FRIConfig = FRIConfig { fasta              :: Text
                           , variants           :: [Variants]
                           , ambiguitycodes     :: [Text]
                           , outputdirectory    :: Text
                           , tsswindowsize      :: Maybe Text
                           , ignorestrandedness :: Bool
                           } deriving (Eq,Show,Read)

data Variants = Variants { vsample     :: Text
                         , vsymbol     :: Text
                         , vchromosome :: Text
                         , vstartpos   :: Text
                         , vendpos     :: Text
                         , vref        :: Text
                         , valt        :: Text
                         , venst       :: Text
                         } deriving (Eq,Show,Read)

instance FromJSON FRIConfig where
  parseJSON (Object v) = parseFRIConfig v
  parseJSON _          = CA.empty 

parseFRIConfig v = FRIConfig
  <$> v .:  "Fasta"
  <*> v .:  "Variants"
  <*> v .:  "Ambiguity_Codes"
  <*> v .:  "Output_Directory"
  <*> v .:? "TSS_Window_Size"
  <*> v .:  "Ignore_Strandness"

instance FromJSON Variants where
  parseJSON (Object v) = parseVariants v
  parseJSON _          = CA.empty

parseVariants v = Variants
  <$> v .: "Sample"
  <*> v .: "HGNC_Symbol"
  <*> v .: "Chromosome"
  <*> v .: "Start_Position"
  <*> v .: "End_Position"
  <*> v .: "Refernce_Allele"
  <*> v .: "Alternate_Allele"
  <*> v .: "ENST"

data BioMartRegion = BioMartRegion { rchromosome :: Text
                                   , rtss        :: Text
                                   , rstrand     :: Text
                                   , rgenename   :: Text
                                   } deriving (Eq,Show,Read)

{--------------------------------------------------------}
