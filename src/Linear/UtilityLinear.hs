{-# LANGUAGE LinearTypes       #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QualifiedDo       #-}

module Linear.UtilityLinear where

import Types
import Utility

import qualified System.IO.Resource.Linear           as Linear
import qualified Control.Functor.Linear              as Control
import           Data.Ord                            as DO
import           Data.Text (Text(..),splitOn,unpack)
import qualified Data.Unrestricted.Linear
import           Data.Vector.Unboxed                 as DVU
import           GHC.Integer
import           Prelude.Linear
import qualified Prelude

getFAILineLinear :: FRIConfig
                 -> BioMartRegion
                 -> Integer
                 -> IO FAI
getFAILineLinear config
                 currentregion
                 seqcounter = do
  cfai <- Linear.run $ Control.do
    fih             <- Linear.openFile (unpack $ fai config)
                                       Linear.ReadMode
    handle'         <- Linear.hSeek fih
                                    Linear.AbsoluteSeek
                                    seqcounter
    (cfai,handle'') <- Linear.hGetLine handle'
    ()              <- Linear.hClose handle''
    Control.return cfai
  case (unpack cfai) == (unpack $ biomartregion_sequencedescription currentregion) of
    True  -> do let cfaif = splitOn "\t"
                                    cfai
                Prelude.return $ FAI { fai_name      = cfaif Prelude.!! 0
                                     , fai_length    = toInteger $
                                                       read      $
                                                       unpack    $
                                                       cfaif Prelude.!! 1
                                     , fai_offset    = toInteger $
                                                       read      $
                                                       unpack    $
                                                       cfaif Prelude.!! 2
                                     , fai_linebases = toInteger $
                                                       read      $
                                                       unpack    $
                                                       cfaif Prelude.!! 3
                                     , fai_linewidth = toInteger $
                                                       read      $
                                                       unpack    $
                                                       cfaif Prelude.!! 4
                                     }
    False -> do let newcharsize = numBytesUtf8String $
                                  unpack cfai
                getFAILineLinear config
                                 currentregion
                                 (plusInteger seqcounter
                                              (toInteger newcharsize)
                                 )

getFASTASequenceLinear :: FRIConfig
                       -> FAI
                       -> Integer
                       -> Vector Char
                       -> IO FASTASequence
getFASTASequenceLinear config
                       cfai
                       seqcounter
                       fastaseq =
  case seqcounter DO.> (fai_length cfai) of
    False -> do newseqchar <- Linear.run $ Control.do
                                ffh            <- Linear.openFile (unpack $ fasta config)
                                                                  Linear.ReadMode
                                handle'        <- Linear.hSeek ffh
                                                               Linear.AbsoluteSeek
                                                               (plusInteger (fai_offset cfai)
                                                                            seqcounter
                                                               )
                                (nsc,handle'') <- Linear.hGetChar handle'
                                ()             <- Linear.hClose handle''
                                Control.return nsc
                let newcharsize = numBytesUtf8Char newseqchar
                case newseqchar of
                  'A' -> getFASTASequenceLinear config 
                                                cfai
                                                (plusInteger seqcounter 
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'T' -> getFASTASequenceLinear config 
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'G' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'C' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'N' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'U' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'R' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'Y' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'K' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'M' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'S' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'W' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'B' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'D' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'H' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  'V' -> getFASTASequenceLinear config
                                                cfai
                                                (plusInteger seqcounter
                                                             (toInteger newcharsize)
                                                )
                                                (DVU.snoc fastaseq
                                                          newseqchar
                                                )
                  _    -> getFASTASequenceLinear config
                                                 cfai
                                                 (plusInteger seqcounter
                                                              (toInteger newcharsize)
                                                 )
                                                 fastaseq
    True  -> Prelude.return $ FASTASequence fastaseq 
