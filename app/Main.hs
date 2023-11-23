{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE Strict            #-}
{-# LANGUAGE TypeOperators     #-}

module Main where

import MacroFont
import RunFRI
import CmdOpts

import Data.Aeson.Types
import Data.List as DL
import Effectful
import Effectful.Ki
import Effectful.Log
import Log.Backend.StandardOutput (withStdOutLogger)
import System.Directory as SD
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO

main :: IO ()
main = do --Print out Fasta Region Inspector ascii art.
          _ <- SIO.putStrLn fri3DMacroFont
          withStdOutLogger $ \logger -> do
            runEff $ runLog "main"
                            logger
                            LogTrace
                            friApp

friApp :: Eff '[Log,IOE] ()
friApp = do --Get command line arguments.
            (args,files) <- liftIO $ SE.getArgs >>= compilerOpts
            --See if files is null.
            if | DL.length files /= 1
               -> do logMessage LogTrace
                                "FRI requires one argument:\n\
                                \Argument 1: Configuration YAML file\n"
                                Null
                     liftIO $ SX.exitWith (SX.ExitFailure 1) 
               | otherwise
               -> do --Run args and files through processArgsAndFiles.
                     _ <- logMessage LogInfo
                                     "Starting up Fasta Region Inspector v0.2.0.0."
                                     Null
                     runStructuredConcurrency $ runFastaRegionInspector (args,files)   
