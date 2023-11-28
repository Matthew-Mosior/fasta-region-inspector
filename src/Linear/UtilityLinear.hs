{-# LANGUAGE LinearTypes       #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QualifiedDo       #-}

module Linear.UtilityLinear where

import Types
import Utility

import qualified System.IO.Resource.Linear           as Linear
import qualified Control.Functor.Linear              as Control
import           Data.List                           as DL
import           Data.List.Linear                    as DLL
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
  case (DLL.takeWhile (/= ' ') (unpack cfai)) == (unpack $ biomartregion_sequencedescription currentregion) of
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
    False -> do newseqline <- Linear.run $ Control.do
                                ffh           <- Linear.openFile (unpack $ fasta config)
                                                                 Linear.ReadMode
                                handle'       <- Linear.hSeek ffh
                                                              Linear.AbsoluteSeek
                                                              (plusInteger (fai_offset cfai)
                                                                           seqcounter
                                                              )
                                (nl,handle'') <- Linear.hGetLine handle'
                                ()            <- Linear.hClose handle''
                                Control.return nl
                case ((toInteger (DVU.length $ fastaseq DVU.++ (DVU.fromList $ unpack newseqline))) DO.> (fai_length cfai)) of
                  True  -> do let newseqlinef  = DL.filter (\x -> x `DL.elem` nucleicacidcodes) $
                                                 unpack newseqline
                              let newseqlineff = DL.take ( (DVU.length $ fastaseq DVU.++ (DVU.fromList $ unpack newseqline))
                                                            -
                                                           (fromIntegral $ fai_length cfai)
                                                         )
                                                 newseqlinef
                              let newlinesize  = numBytesUtf8String newseqlineff
                              getFASTASequenceLinear config
                                                     cfai
                                                     ( plusInteger seqcounter
                                                                   (toInteger newlinesize)
                                                     )
                                                     ( fastaseq
                                                       DVU.++
                                                       (DVU.fromList newseqlinef)
                                                     )
                  False -> do let newseqlinef  = DL.filter (\x -> x `DL.elem` nucleicacidcodes) $
                                                 unpack newseqline
                              let newlinesize  = numBytesUtf8String newseqlinef
                              getFASTASequenceLinear config
                                                     cfai
                                                     ( plusInteger seqcounter
                                                                   (toInteger newlinesize)
                                                     )
                                                     ( fastaseq
                                                       DVU.++
                                                       (DVU.fromList newseqlinef)
                                                     ) 
    True  -> Prelude.return $ FASTASequence fastaseq
  where
    nucleicacidcodes = [ 'A'
                       , 'T'
                       , 'G'
                       , 'C'
                       , 'N'
                       , 'U'
                       , 'R'
                       , 'Y'
                       , 'K'
                       , 'M'
                       , 'S'
                       , 'W'
                       , 'B'
                       , 'D'
                       , 'H'
                       , 'V'
                       ]
