{-=Fasta-Region-Inspector (FRI): A Somatic=-}
{-=Hypermutation Analysis Tool.=-}
{-=Author: Matthew Mosior=-}
{-=Synposis: This Haskell script will=-}
{-=process command line arguments to FRI.=-}


{-Language extension.-}

{-# LANGUAGE Strict            #-}
{-# LANGUAGE StrictData        #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE OverloadedLists   #-} 

{---------------------}


{-Module.-}

module QueryBioMart where

{---------}


{-Import modules.-}

import YamlParser

{-----------------}


{-Imports.-}

import Data.List as DL
import Data.List.Split as DLS
import Data.Text as DText
import qualified Data.Text.Lazy as DTL
import qualified Data.Text.Lazy.IO as DTLIO
import qualified Data.Map as Map
import Data.Time as DTime
import System.Console.GetOpt as SCG
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO
import System.Process as SP
import Text.XML

{----------}


{-Query BioMart functions.-}

runQueryBioMart :: FRIConfig -> IO [BioMartRegion]
runQueryBioMart config = do
  --Generate BioMart compatible xml.
  currenttandd <- DTime.getZonedTime 
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Generating BioMart compatible XML ...")
  
  let biomartxml = Element "Query" [ ("virtualSchemaName","Default")
                                   , ("formatter","CSV")
                                   , ("header","0")
                                   , ("uniqueRows","0")
                                   , ("count","")
                                   , ("datasetConfigVersion","0.6")
                                   ]
                     [ NodeElement $ Element "Dataset" [ ("name","hsapiens_gene_ensembl")
                                                       , ("interface","default")
                                                       ]
                       [ NodeElement $ Element "Filter" [ ("name","ensembl_transcript_id")
                                                        , ("value",finalallensts)
                                                        ] [],
                         NodeElement $ Element "Attribute" [ ("name","chromosome_name")
                                                           ] [] ,
                         NodeElement $ Element "Attribute" [ ("name","transcription_start_site")
                                                           ] [] ,
                         NodeElement $ Element "Attribute" [ ("name","strand")
                                                           ] [] ,
                         NodeElement $ Element "Attribute" [ ("name","external_gene_name")
                                                           ] []
                       ]
                     ]
  let finalbiomartxml = Document (Prologue []
                                           (Just Doctype { doctypeName = "Query" , doctypeID = Nothing })
                                           [])
                                 biomartxml
                                 []
  --Query BioMart using wget command via system process call.
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                         ++ "Querying and downloading region data from BioMart via system process call ...")
  (_,hout,_,ph) <- SP.createProcess $
                     (SP.shell ("wget " ++
                                "'http://www.ensembl.org/biomart/martservice?query=" ++
                                DText.unpack (DTL.toStrict (renderText def $ finalbiomartxml))
                     )) -- { std_out = CreatePipe }
  ec <- waitForProcess ph
  case ec of
    SX.ExitFailure _ -> do currenttandd <- DTime.getZonedTime
                           _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                                                  ++ "Could not query and download region data from BioMart via system process call.")
                           return []
    SX.ExitSuccess   -> do currenttandd <- DTime.getZonedTime
                           _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                                                  ++ "Successfully queried and downloaded region data from BioMart via system process call ...")
                           --Process downloaded BioMart data.
                           currenttandd <- DTime.getZonedTime
                           _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                                                  ++ "Beginning to process downloaded region data ...") 
                           case hout of
                             Just res -> do --Turn handle into usable [BioMartRegion].
                                            currenttandd <- DTime.getZonedTime
                                            _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                                                                   ++ "Converting std out handle from query to region data ...")
                                            returnedbiomartregions <- SIO.hGetContents res 
                                            let finalbiomartregions    = DL.map (\x -> DLS.splitOn "," x)
                                                                         (DL.lines returnedbiomartregions)
                                            let processedregiondata    = DL.concat
                                                                         (DL.map 
                                                                         (DL.map (\y -> BioMartRegion { rchromosome = DText.pack
                                                                                                                     (DL.take 1 y)
                                                                                                     , rtss        = DText.pack
                                                                                                                     (DL.drop 1
                                                                                                                     (DL.take 2 y))
                                                                                                     , rstrand     = DText.pack
                                                                                                                     (DL.drop 2
                                                                                                                     (DL.take 3 y))
                                                                                                     , rgenename   = DText.pack
                                                                                                                     (DL.drop 3
                                                                                                                     (DL.take 4 y))
                                                                                                     }))
                                                                         finalbiomartregions)
                                            return processedregiondata
                             Nothing  -> do --Unable to turn std out handl into usable [BioMartRegion].
                                            currenttandd <- DTime.getZonedTime
                                            _ <- SIO.putStrLn ("[" ++ (show currenttandd) ++ "] "
                                                                   ++ "Unable to convert std out handle from query to region data ...")
                                            return []
    where
      finalallensts = DText.pack
                      (DL.intercalate
                      ","
                      allensts)
      allensts      = DL.nub
                      (DL.map (\x -> DText.unpack x)
                      (DL.map (\x -> venst x)
                      allvariants))
      allvariants   = variants config

{--------------------------}
