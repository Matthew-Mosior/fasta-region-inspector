{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE Strict            #-}

module Main where

import CmdOpts
import Logging
import MacroFont
import RunFRI

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
               -> liftIO $ SX.exitWith (SX.ExitFailure 1) 
               | otherwise
               -> runFRILogging $ runStructuredConcurrency $ runFastaRegionInspector (args,files)   
