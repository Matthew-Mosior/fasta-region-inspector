{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE StrictData    #-}

module Types where

import Data.ByteString
import Data.List as DL
import Data.Maybe as DMaybe
import Data.Text as DText
import Control.Applicative as CA
import GHC.Generics 

data FRIConfig = FRIConfig { fasta                :: Text
                           , fai                  :: Text
                           , variants             :: [Variant]
                           , ambiguitycodes       :: [Text]
                           , outputdirectory      :: Text
                           , tsswindowsize        :: Maybe Text
                           , keepbiomart          :: Bool
                           , ignorestrandedness   :: Bool
                           , writeambiguitycodes  :: Bool
                           , maxnumberconcthreads :: Int
                           } deriving (Eq,Generic,Read,Show)

data Variant = Variant { variant_sample              :: Text
                       , variant_symbol              :: Text
                       , variant_sequencedescription :: Text
                       , variant_startpos            :: Text
                       , variant_endpos              :: Text
                       , variant_ref                 :: Text
                       , variant_alt                 :: Text
                       , variant_enst                :: Text
                       } deriving (Eq,Generic,Read,Show)

data BioMartRegion = BioMartRegion { biomartregion_sequencedescription :: Text
                                   , biomartregion_tss                 :: Text
                                   , biomartregion_strand              :: Text
                                   , biomartregion_genename            :: Text
                                   } deriving (Eq,Generic,Read,Show)

data FAI = FAI { fai_name      :: Text
               , fai_length    :: Integer
               , fai_offset    :: Integer
               , fai_linebases :: Integer
               , fai_linewidth :: Integer
               } deriving (Eq,Generic,Read,Show)

newtype FASTASequence = FASTASequence ByteString
  deriving (Eq,Generic,Read,Show)
