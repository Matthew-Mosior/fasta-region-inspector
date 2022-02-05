{-=Fasta-Region-Inspector (FRI): A Somatic=-}
{-=Hypermutation Analysis Tool.=-}
{-=Author: Matthew Mosior=-}
{-=Synposis: This Haskell script will=-}
{-=launch the FRI application.=-}

{-Language extension.-}

{-# LANGUAGE Strict     #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE MultiWayIf #-}

{---------------------}


{-Module.-}

module Main where

{---------}


{-Import src.-}

import RunFRI
import CmdOpts

{-------------}


{-Imports.-}

import Data.List as DL
import System.Directory as SD
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO

{----------}


{-Main function.-}

main :: IO ()
main = do
  --Get command line arguments.
  (args,files) <- SE.getArgs >>= compilerOpts
  --See if files is null.
  if | (DL.length files) /= 1
     -> do --Print error statement and exit.
           SIO.putStrLn "FRI requires one argument:\n\
                        \Argument 1: Configuration YAML file\n"
           SX.exitWith (SX.ExitFailure 1)
     | otherwise
     -> do --Run args and files through processArgsAndFiles.
           runFastaRegionInspector (args,files)

{----------------}
