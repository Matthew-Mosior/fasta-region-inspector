{-# LANGUAGE LinearTypes       #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QualifiedDo       #-}

module Linear.UtilityLinear where

import Types
import Utility

import qualified System.IO.Resource.Linear           as Linear
import qualified Control.Functor.Linear              as Control
import qualified Data.ByteString.Char8               as DB
import           Data.List                           as DL
import           Data.List.Linear                    as DLL
import           Data.Ord                            as DO
import           Data.Text (Text(..),splitOn,unpack)
import qualified Data.Unrestricted.Linear
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

getFASTASequenceLinearReverse :: FRIConfig
                              -> BioMartRegion
                              -> FAI
                              -> DB.ByteString
                              -> Linear.Handle %1
                              -> Linear.RIO (Ur (Maybe FASTASequence))
getFASTASequenceLinearReverse config
                              currentregion
                              cfai
                              fastaseq
                              ffh =
  loop config
       currentregion
       cfai
       fastaseq
       ffh
    where
      loop :: FRIConfig
           -> BioMartRegion
           -> FAI
           -> DB.ByteString
           -> Linear.Handle %1
           -> Linear.RIO (Ur (Maybe FASTASequence))
      loop config
           currentregion
           cfai
           fastaseq
           ffh = Control.do
        (Ur isEOF,handle') <- Linear.hIsEOF ffh
        case isEOF of
          True  -> Control.do () <- Linear.hClose handle'
                              Control.return $ Ur Nothing
          False -> Control.do (Ur newseqline,handle'') <- Linear.hGetLine handle'
                              let newseqlinef = DB.append (DB.pack $ unpack newseqline)
                                                          fastaseq
                              case (tsswindowsize config) of
                                Just size -> case (DB.length newseqlinef) DO.> (read $ unpack size) of
                                               True  -> Control.do let newseqlineff  = DB.dropEnd ( (DB.length newseqlinef)
                                                                                                    -
                                                                                                    (read $ unpack size)
                                                                                                  )
                                                                                       newseqlinef
                                                                   ()               <- Linear.hClose handle''
                                                                   Control.return $ Ur $ Just $ FASTASequence newseqlineff
                                               False -> Control.do loop config
                                                                        currentregion
                                                                        cfai
                                                                        newseqlinef
                                                                        handle''
                                Nothing   -> case (DB.length newseqlinef) DO.> 2000 of
                                               True  -> Control.do let newseqlineff  = DB.dropEnd ( (DB.length newseqlinef)
                                                                                                    -
                                                                                                    2000
                                                                                                  )
                                                                                       newseqlinef
                                                                   ()               <- Linear.hClose handle''
                                                                   Control.return $ Ur $ Just $ FASTASequence newseqlineff
                                               False -> Control.do loop config
                                                                        currentregion
                                                                        cfai
                                                                        newseqlinef
                                                                        handle''

getFASTASequenceLinearForward :: FRIConfig
                              -> BioMartRegion
                              -> FAI
                              -> DB.ByteString
                              -> Linear.Handle %1
                              -> Linear.RIO (Ur (Maybe FASTASequence))
getFASTASequenceLinearForward config
                              currentregion
                              cfai
                              fastaseq
                              ffh =
  loop config
       currentregion
       cfai
       fastaseq
       ffh
    where
      loop :: FRIConfig
           -> BioMartRegion
           -> FAI
           -> DB.ByteString
           -> Linear.Handle %1
           -> Linear.RIO (Ur (Maybe FASTASequence))
      loop config
           currentregion
           cfai
           fastaseq
           ffh = Control.do
        (Ur isEOF,handle') <- Linear.hIsEOF ffh
        case isEOF of
          True  -> Control.do () <- Linear.hClose handle'
                              Control.return $ Ur Nothing
          False -> Control.do (Ur newseqline,handle'') <- Linear.hGetLine handle'
                              let newseqlinef = DB.append fastaseq
                                                          (DB.pack $ unpack newseqline)
                              case (tsswindowsize config) of
                                Just size -> case (DB.length newseqlinef) DO.> (read $ unpack size) of
                                               True  -> Control.do let newseqlineff  = DB.dropEnd ( (DB.length newseqlinef)
                                                                                                    -
                                                                                                    (read $ unpack size)
                                                                                                  )
                                                                                       newseqlinef
                                                                   ()               <- Linear.hClose handle''
                                                                   Control.return $ Ur $ Just $ FASTASequence newseqlineff
                                               False -> Control.do loop config
                                                                        currentregion
                                                                        cfai
                                                                        newseqlinef
                                                                        handle''
                                Nothing   -> case (DB.length newseqlinef) DO.> 2000 of
                                               True  -> Control.do let newseqlineff  = DB.dropEnd ( (DB.length newseqlinef)
                                                                                                    -
                                                                                                    2000
                                                                                                  )
                                                                                       newseqlinef
                                                                   ()               <- Linear.hClose handle''
                                                                   Control.return $ Ur $ Just $ FASTASequence newseqlineff
                                               False -> Control.do loop config
                                                                        currentregion
                                                                        cfai
                                                                        newseqlinef
                                                                        handle''

getFASTASequenceLinear :: FRIConfig
                       -> BioMartRegion
                       -> FAI
                       -> Linear.RIO (Ur (Maybe FASTASequence))
getFASTASequenceLinear config
                       currentregion
                       cfai = 
  if | currentregionstrand == "-1"
     -> Control.do handle  <- Linear.openFile (unpack $ fasta config)
                                              Linear.ReadMode
                   handle' <- case (tsswindowsize config) of
                                Just size -> Linear.hSeek handle
                                                          Linear.AbsoluteSeek
                                                          ( minusInteger ( plusInteger (fai_offset cfai)
                                                                                       ( plusInteger (toInteger $ read $ unpack $ biomartregion_tss currentregion)
                                                                                                     1
                                                                                       )
                                                                         )
                                                                         (read $ unpack size)
                                                          )
                                Nothing   -> Linear.hSeek handle
                                                          Linear.AbsoluteSeek
                                                          ( minusInteger ( plusInteger (fai_offset cfai)
                                                                                       ( plusInteger (toInteger $ read $ unpack $ biomartregion_tss currentregion)
                                                                                                     1
                                                                                       )
                                                                         ) 2000
                                                          )
                   getFASTASequenceLinearReverse config
                                                 currentregion
                                                 cfai
                                                 DB.empty
                                                 handle'
     | otherwise
     -> Control.do handle  <- Linear.openFile (unpack $ fasta config)
                                              Linear.ReadMode
                   handle' <- Linear.hSeek handle
                                            Linear.AbsoluteSeek
                                            ( plusInteger (fai_offset cfai) ( plusInteger (toInteger $ read $ unpack $ biomartregion_tss currentregion)
                                                                                          1
                                                                            )
                                            )
                   getFASTASequenceLinearForward config
                                                 currentregion
                                                 cfai
                                                 DB.empty
                                                 handle'
    where
      currentregionstrand = unpack $ biomartregion_strand currentregion
