{-=Fasta-Region-Inspector (FRI): A Somatic=-}
{-=Hypermutation Analysis Tool.=-}
{-=Author: Matthew Mosior=-}
{-=Synposis: This Haskell script will=-}
{-=process command line arguments to FRI.=-}


{-Module.-}

module CmdOpts where

{---------}


{-Imports.-}

import Data.List as DL
import System.Console.GetOpt as SCG
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO

{----------}


{-Custom CML Option Datatype.-}

data Flag = Help -- --help
  deriving (Eq,Ord,Show)

{-----------------------------}


{-Option datatype function relating to datatype above.-}

--options -> This function will
--describe flags.
options :: [OptDescr Flag]
options = 
  [ Option ['h'] ["help"] (NoArg Help) "Print this help message.\n"
  ]

{------------------------------------------------------}


{-Function to correctly parse the flags.-}

compilerOpts :: [String] -> IO ([Flag],[String])
compilerOpts argv =
  case getOpt Permute options argv of
    (args,files,[]) ->
      if | DL.elem Help args
         -> do SIO.hPutStrLn stderr (greeting ++ SCG.usageInfo header options)
               SX.exitWith SX.ExitSuccess
         | otherwise
         -> return (DL.nub args,files)
    (_,_,errors) -> do SIO.hPutStrLn stderr (DL.concat errors ++ SCG.usageInfo header options)
                       SX.exitWith (SX.ExitFailure 1)
    where
      greeting = "Fasta-Region-Inspector\n\
                 \Copyright (c) Matthew C. Mosior 2022\n"
      header   = "Usage: stack exec fasta-region-inspector-exe [-h] [Configuration YAML]\n\
                 \Fasta-Region-Inspector (FRI), Version 0.1.0.0\n"
      
{----------------------------------------}
