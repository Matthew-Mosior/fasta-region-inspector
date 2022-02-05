{-=Fasta-Region-Inspector (FRI): A Somatic=-}
{-=Hypermutation Analysis Tool.=-}
{-=Author: Matthew Mosior=-}
{-=Synposis: This Haskell script will=-}
{-=launch the FRI main runner function.=-}


{-Language extension.-}

{-# LANGUAGE Strict      #-}
{-# LANGUAGE StrictData  #-}
{-# LANGUAGE MultiWayIf  #-}
{-# LANGUAGE QuasiQuotes #-}

{---------------------}


{-Module.-}

module RunFRI where

{---------}


{-Import modules.-}

import CmdOpts
import Common
import Inspect
import YamlParser

{-----------------}


{-Imports.-}

import Data.ByteString.Char8 as DBC
import Data.Char as DC
import Data.List as DL
import Data.Maybe as DMaybe
import Data.Text as DText
import Data.Yaml as DYaml
import System.Directory as SD
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO
import Text.RawString.QQ as TRQQ
import Text.Regex as TR
import Text.Regex.TDFA as TRTDFA

{----------}


{-3D Macro Font.-}

fri3DMacroFont :: String
fri3DMacroFont =
  [r|
        ______           __           ____             _                ____                           __            
       / ____/___ ______/ /_____ _   / __ \___  ____ _(_)___  ____     /  _/___  _________  ___  _____/ /_____  _____
      / /_  / __ `/ ___/ __/ __ `/  / /_/ / _ \/ __ `/ / __ \/ __ \    / // __ \/ ___/ __ \/ _ \/ ___/ __/ __ \/ ___/
     / __/ / /_/ (__  ) /_/ /_/ /  / _, _/  __/ /_/ / / /_/ / / / /  _/ // / / (__  ) /_/ /  __/ /__/ /_/ /_/ / /    
    /_/    \__,_/____/\__/\__,_/  /_/ |_|\___/\__, /_/\____/_/_/_/__/___/_/ /_/____/ .___/\___/\___/\__/\____/_/     
                                             /____/_/ __ \ <  // __ \ / __ \      /_/                                
                                             | | / / / / / / // / / // / / /                                         
                                             | |/ / /_/ / / // /_/ // /_/ /                                          
                                             |___/\____(_)_(_)____(_)____/                                           
                                                                  
                                           Copyright (c) Matthew C. Mosior 2022          
  |]

{----------------}


{-FRI specific functions.-}

tssWindowSizeCheck :: FRIConfig -> Bool
tssWindowSizeCheck config = if | DMaybe.isJust (tsswindowsize config)
                               -> if | DL.all (DC.isDigit)
                                              (DText.unpack
                                              (DMaybe.fromJust
                                              (tsswindowsize config)))
                                     -> True
                                     | otherwise
                                     -> False
                               | otherwise
                               -> True

ambiguityCodesCheck :: FRIConfig -> Bool
ambiguityCodesCheck config = if | DL.all 
                                  (DL.all (\y -> y `DL.elem` nucleotideambiguitycodes)) 
                                          (DL.map (DText.unpack)
                                          (ambiguitycodes config))
                                -> True
                                | otherwise
                                -> False 
  where
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

processConfigurationYaml :: FRIConfig -> IO Bool
processConfigurationYaml config = do
  doesinputfastafileexist <- SD.doesFileExist (DText.unpack $ fasta config)
  if | doesinputfastafileexist &&
       tssWindowSizeCheck config &&
       ambiguityCodesCheck config
     -> return True
     | otherwise
     -> return False  


runFastaRegionInspector :: ([Flag],[String]) -> IO ()
runFastaRegionInspector ([],[]) = return ()
runFastaRegionInspector (options,inputfiles) = do
  --Print out Fasta Region Inspector ascii art.
  SIO.putStrLn fri3DMacroFont
  --Read in configuration YAML.
  readinputyaml <- DBC.readFile ((\(x:_) -> x) inputfiles)
  --Decode readinputyaml.
  decodedinputyaml <- decodeThrow readinputyaml :: IO FRIConfig
  --Process and ensure correctly formatting Configuration YAML input.
  processedConfig <- processConfigurationYaml decodedinputyaml
  if | processedConfig
     -> do --Run fastaRegionInspect decodedinputyaml
           fastaRegionInspect decodedinputyaml
     | otherwise
     -> do --Print out failure message.
           SIO.putStrLn "Could not sanitize Configuration YAML.\n"
           SX.exitWith (SX.ExitFailure 1) 

{-------------------------}
