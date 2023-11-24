{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE Strict            #-}
{-# LANGUAGE TypeOperators     #-}

module Main where

import Logging
import MacroFont
import RunFRI
import CmdOpts

import Data.Aeson.Types
import Data.List as DL
import Effectful
import Effectful.Ki
import System.Directory as SD
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO

main :: IO ()
main = do --Print out Fasta Region Inspector ascii art.
          _ <- SIO.putStrLn fri3DMacroFont
          runEff friApp

friApp :: Eff '[IOE] ()
friApp = do --Get command line arguments.
            (args,files) <- liftIO $ SE.getArgs >>= compilerOpts
            --See if files is null.
            if | DL.length files /= 1
               -> do _ <- liftIO $ showPrettyLog LogTrace
                                                 "friApp"
                                                 "FRI requires one argument: Configuration YAML file."
                     liftIO $ SX.exitWith (SX.ExitFailure 1) 
               | otherwise
               -> do _ <- liftIO $ showPrettyLog LogInfo
                                                 "friApp"
                                                 "Starting up fasta-region-inspector v0.2.0.0." 
                     runStructuredConcurrency $ runFastaRegionInspector (args,files)   
