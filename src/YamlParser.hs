{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE OverloadedStrings #-}

module YamlParser where

import Types

import Data.Aeson   as DAeson
import Data.HashMap.Lazy as DHL
import Data.List as DL
import Data.Maybe as DMaybe
import Data.Text as DText
import Data.Yaml as DYaml
import Control.Applicative as CA
import GHC.Generics 

instance FromJSON FRIConfig where
  parseJSON (Object v) = parseFRIConfig v
  parseJSON _          = CA.empty 

parseFRIConfig v = FRIConfig
  <$> v .:  "Fasta"
  <*> v .:  "Fasta_Index"
  <*> v .:  "Variants"
  <*> v .:  "Ambiguity_Codes"
  <*> v .:  "Output_Directory"
  <*> v .:? "TSS_Window_Size"
  <*> v .:  "Keep_BioMart"
  <*> v .:  "Ignore_Strandedness"
  <*> v .:  "Write_Ambiguity_Codes"

instance FromJSON Variant where
  parseJSON (Object v) = parseVariant v
  parseJSON _          = CA.empty

parseVariant v = Variant
  <$> v .: "Sample"
  <*> v .: "HGNC_Symbol"
  <*> v .: "Sequence_Description"
  <*> v .: "Start_Position"
  <*> v .: "End_Position"
  <*> v .: "Reference_Allele"
  <*> v .: "Alternate_Allele"
  <*> v .: "ENST"
