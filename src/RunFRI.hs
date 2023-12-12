{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE KindSignatures    #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE RankNTypes        #-}
{-# LANGUAGE Strict            #-}
{-# LANGUAGE TypeOperators     #-}

module RunFRI where

import CmdOpts
import Inspect
import Logging
import Types
import YamlParser

import Data.ByteString.Char8 as DBC
import Data.Char as DC
import Data.List as DL
import Data.Maybe as DMaybe
import Data.Text as DText
import Data.Yaml as DYaml
import Effectful
import Effectful.Ki
import GHC.Stack                    (callStack,getCallStack,HasCallStack)
import System.Directory as SD
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO
import Text.Regex as TR
import Text.Regex.TDFA as TRTDFA

tssWindowSizeCheck :: FRIConfig
                   -> Bool
tssWindowSizeCheck config =
  if | DMaybe.isJust (tsswindowsize config)
     -> if | DL.all (DC.isDigit)
                    (DText.unpack
                    (DMaybe.fromJust
                    (tsswindowsize config)))
           -> True
           | otherwise
           -> False
     | otherwise
     -> True

ambiguityCodesCheck :: FRIConfig
                    -> Bool
ambiguityCodesCheck config =
  if | DL.all 
       (DL.all (\y -> y `DL.elem` nucleotideambiguitycodes)) 
               (DL.map (DText.unpack)
               (ambiguitycodes config))
     -> True
     | otherwise
     -> False 
  where
    nucleotideambiguitycodes :: [Char]
    nucleotideambiguitycodes = ['A'
                               ,'G'
                               ,'C'
                               ,'T'
                               ,'Y'
                               ,'R'
                               ,'W'
                               ,'S'
                               ,'K'
                               ,'M'
                               ,'D'
                               ,'V'
                               ,'H'
                               ,'B'
                               ,'X'
                               ,'N'
                               ,'-']

maxNumberofConcurrentThreadsSizeCheck :: FRIConfig
                                      -> Bool
maxNumberofConcurrentThreadsSizeCheck config = 
  if | maxnumberconcthreads config < 1000 
     -> True
     | otherwise
     -> False

processConfigurationYaml :: FRIConfig
                         -> IO Bool
processConfigurationYaml config = do
  doesinputfastafileexist  <- SD.doesFileExist (DText.unpack $ fasta config)
  doesinputfastaindexexist <- SD.doesFileExist (DText.unpack $ fai config)
  if | doesinputfastafileexist &&
       doesinputfastaindexexist &&
       tssWindowSizeCheck config &&
       ambiguityCodesCheck config &&
       maxNumberofConcurrentThreadsSizeCheck config
     -> return True
     | otherwise
     -> return False  

runFastaRegionInspector :: forall {es :: [Effect]} {b}.
                           ( FRILogging :> es
                           , StructuredConcurrency :> es
                           , IOE :> es
                           , HasCallStack
                           )
                        => ([Flag],[FilePath])
                        -> Eff es ()
runFastaRegionInspector ([],[])        = liftIO $ return ()
runFastaRegionInspector (_,inputfiles) = do
  --Get calling function name.
  let ((callingfunction,_):_) = getCallStack callStack
  --Read in configuration YAML.
  readinputyaml <- liftIO $ DBC.readFile ((\(x:_) -> x) inputfiles)
  --Decode readinputyaml.
  decodedinputyaml <- decodeThrow readinputyaml -- :: IO FRIConfig
  --Process and ensure correctly formatting Configuration YAML input.
  processedConfig <- liftIO $ processConfigurationYaml decodedinputyaml
  if | processedConfig
        --Process the input yaml.
     -> fastaRegionInspect decodedinputyaml
     | otherwise
     -> do --Print out failure message.
           _ <- showPrettyLog LogDebug
                              (maxnumberconcthreads decodedinputyaml)
                              callingfunction
                              "Could not sanitize Configuration YAML."
           liftIO $ SX.exitWith (SX.ExitFailure 1)
