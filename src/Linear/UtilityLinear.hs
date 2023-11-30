{-# LANGUAGE LinearTypes       #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}
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

getFAILineLinearS :: FRIConfig
                  -> BioMartRegion
                  -> Linear.Handle %1
                  -> Linear.RIO (Ur (Maybe FAI))
getFAILineLinearS config
                  currentregion
                  fih =
  loop config
       currentregion
       fih
    where
      loop :: FRIConfig
           -> BioMartRegion
           -> Linear.Handle %1
           -> Linear.RIO (Ur (Maybe FAI))
      loop config
           currentregion
           fih = Control.do
        (Ur isEOF,handle') <- Linear.hIsEOF fih
        case isEOF of
          True  -> Control.do () <- Linear.hClose handle'
                              Control.return $ Ur Nothing
          False -> Control.do (Ur cfai,handle'') <- Linear.hGetLine handle'
                              case (DLL.takeWhile (/= '\t') (unpack cfai)) == (unpack $ biomartregion_sequencedescription currentregion) of
                                True  -> Control.do () <- Linear.hClose handle''
                                                    let cfaif = splitOn "\t"
                                                                        cfai
                                                    Control.return $ Ur $ Just FAI { fai_name      = cfaif Prelude.!! 0
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
                                False -> Control.do loop config
                                                         currentregion
                                                         handle''

getFAILineLinear :: FRIConfig
                 -> BioMartRegion
                 -> Linear.RIO (Ur (Maybe FAI))
getFAILineLinear config
                 currentregion = Control.do
  handle <- Linear.openBinaryFile (unpack $ fai config)
                                  Linear.ReadMode
  getFAILineLinearS config
                    currentregion
                    handle

getFASTASequenceLinearS :: FRIConfig
                        -> FAI
                        -> Integer
                        -> Vector Char
                        -> Linear.Handle %1
                        -> Linear.RIO (Ur (Maybe FASTASequence))
getFASTASequenceLinearS config
                        cfai
                        seqcounter
                        fastaseq 
                        ffh =
  loop config
       cfai
       seqcounter
       fastaseq
       ffh
    where
      loop :: FRIConfig
           -> FAI
           -> Integer
           -> Vector Char
           -> Linear.Handle %1
           -> Linear.RIO (Ur (Maybe FASTASequence))
      loop config
           cfai
           seqcounter
           fastaseq
           ffh = Control.do
        (Ur isEOF,handle') <- Linear.hIsEOF ffh
        case isEOF of
          True  -> Control.do () <- Linear.hClose handle'
                              Control.return $ Ur Nothing
          False -> Control.do (Ur newseqline,handle'') <- Linear.hGetLine handle'
                              case seqcounter DO.> (fai_length cfai) of
                                False -> case ((toInteger (DVU.length $ fastaseq DVU.++ (DVU.fromList $ unpack newseqline))) DO.> (fai_length cfai)) of
                                           True  -> Control.do let newseqlinef  = DL.take ( (DVU.length $ fastaseq DVU.++ (DVU.fromList $ unpack newseqline))
                                                                                             -
                                                                                            (fromIntegral $ fai_length cfai)
                                                                                          ) $
                                                                                   unpack newseqline
                                                               let newlinesize  = numBytesUtf8String newseqlinef
                                                               ()               <- Linear.hClose handle''
                                                               Control.return $ Ur $ Just $ FASTASequence $ DVU.fromList newseqlinef
                                           False -> Control.do let newlinesize  = numBytesUtf8String $
                                                                                  unpack newseqline 
                                                               loop config
                                                                    cfai
                                                                    ( plusInteger seqcounter
                                                                                  (toInteger newlinesize)
                                                                    )
                                                                    ( fastaseq
                                                                      DVU.++
                                                                      (DVU.fromList $ unpack newseqline)
                                                                    )
                                                                    handle''
                                True  -> Control.do () <- Linear.hClose handle''
                                                    Control.return $ Ur $ Just $ FASTASequence fastaseq

getFASTASequenceLinear :: FRIConfig
                       -> FAI
                       -> Linear.RIO (Ur (Maybe FASTASequence))
getFASTASequenceLinear config
                       cfai = Control.do
  handle <- Linear.openBinaryFile (unpack $ fasta config)
                                  Linear.ReadMode
  getFASTASequenceLinearS config
                          cfai
                          0
                          DVU.empty
                          handle
