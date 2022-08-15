{-=Fasta-Region-Inspector (FRI): A Somatic=-}
{-=Hypermutation Analysis Tool.=-}
{-=Author: Matthew Mosior=-}
{-=Synposis: This Haskell script will=-}
{-=process command line arguments to FRI.=-}


{-Module.-}

module QueryBioMart where

{---------}


{-Import modules.-}

import Common
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
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Generating BioMart compatible XML ...")
  
  let biomartxml = Element "Query" [ ("virtualSchemaName","default")
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
  _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                         ++ "Querying and downloading region data from BioMart via system process call ...") 
  if | DL.last outputdir == '/'
     -> do (_,_,_,ph) <- SP.createProcess $
                              (SP.shell ("wget "                                                        ++
                                         "-O "                                                          ++
                                         outputdir                                                      ++
                                         "biomartresult.txt "                                           ++
                                         "'http://www.ensembl.org/biomart/martservice?query="           ++
                                         DText.unpack (DTL.toStrict (renderText def $ finalbiomartxml)) ++
                                         "'"
                              )) 
           ec <- waitForProcess ph
           case ec of
             SX.ExitFailure _ -> do currenttandd <- DTime.getZonedTime
                                    _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                           ++ "Could not query and download region data from BioMart via system process call.")
                                    return []
             SX.ExitSuccess   -> do currenttandd <- DTime.getZonedTime
                                    _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                           ++ "Successfully queried and downloaded region data from BioMart via system process call ...") 
                                    currenttandd <- DTime.getZonedTime
                                    _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                           ++ "Beginning to process downloaded region data ...")
                                    returnedbiomartregions <- SIO.readFile' (outputdir ++
                                                                             "biomartresult.txt")
                                    let finalbiomartregions = DL.map (\x -> DLS.splitOn "," x)
                                                                     (DL.lines returnedbiomartregions) 
                                    let processedregiondata    = DL.map (\x -> BioMartRegion { rchromosome = DText.pack
                                                                                                             (DL.concat
                                                                                                             (DL.take 1 x))
                                                                                             , rtss        = DText.pack
                                                                                                             (DL.concat
                                                                                                             (DL.drop 1
                                                                                                             (DL.take 2 x)))
                                                                                             , rstrand     = DText.pack
                                                                                                             (DL.concat
                                                                                                             (DL.drop 2
                                                                                                             (DL.take 3 x)))
                                                                                             , rgenename   = DText.pack
                                                                                                             (DL.concat
                                                                                                             (DL.drop 3
                                                                                                             (DL.take 4 x)))
                                                                                             }
                                                                        )
                                                                 finalbiomartregions
                                    if | keepbiomart config
                                       -> return processedregiondata
                                       | otherwise
                                       -> do currenttandd <- DTime.getZonedTime
                                             _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                                    ++ "Removing downloaded BioMart region file ...")
                                             (_,_,_,ph) <- SP.createProcess $
                                                                (SP.shell ("rm "     ++
                                                                           outputdir ++
                                                                           "biomartresult.txt"
                                                                           ))
                                             ec <- waitForProcess ph
                                             case ec of
                                               SX.ExitFailure _ -> do currenttandd <- DTime.getZonedTime
                                                                      _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                                                             ++ "Could not remove downloaded BioMart region file ...")
                                                                      return processedregiondata
                                               SX.ExitSuccess   -> do currenttandd <- DTime.getZonedTime
                                                                      _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                                                             ++ "Removed downloaded BioMart region file ...")
                                                                      return processedregiondata
     | otherwise
     -> do (_,_,_,ph) <- SP.createProcess $
                              (SP.shell ("wget "                                                        ++
                                         "-O "                                                          ++
                                         outputdir                                                      ++
                                         "/"                                                            ++
                                         "biomartresult.txt "                                           ++
                                         "'http://www.ensembl.org/biomart/martservice?query="           ++
                                         DText.unpack (DTL.toStrict (renderText def $ finalbiomartxml)) ++
                                         "'"
                              )) 
           ec <- waitForProcess ph
           case ec of
             SX.ExitFailure _ -> do currenttandd <- DTime.getZonedTime
                                    _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                           ++ "Could not query and download region data from BioMart via system process call.")
                                    return []
             SX.ExitSuccess   -> do currenttandd <- DTime.getZonedTime
                                    _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                           ++ "Successfully queried and downloaded region data from BioMart via system process call ...") 
                                    currenttandd <- DTime.getZonedTime
                                    _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                           ++ "Beginning to process downloaded region data ...")
                                    returnedbiomartregions <- SIO.readFile' (outputdir ++
                                                                             "/"       ++
                                                                             "biomartresult.txt")
                                    let finalbiomartregions = DL.map (\x -> DLS.splitOn "," x)
                                                                     (DL.lines returnedbiomartregions)                                   
                                    let processedregiondata    = DL.map (\x -> BioMartRegion { rchromosome = DText.pack
                                                                                                             (DL.concat
                                                                                                             (DL.take 1 x))
                                                                                             , rtss        = DText.pack
                                                                                                             (DL.concat
                                                                                                             (DL.drop 1
                                                                                                             (DL.take 2 x)))
                                                                                             , rstrand     = DText.pack
                                                                                                             (DL.concat
                                                                                                             (DL.drop 2
                                                                                                             (DL.take 3 x)))
                                                                                             , rgenename   = DText.pack
                                                                                                             (DL.concat
                                                                                                             (DL.drop 3
                                                                                                             (DL.take 4 x)))
                                                                                             }
                                                                        )
                                                                 finalbiomartregions  
                                    if | keepbiomart config
                                       -> return processedregiondata
                                       | otherwise
                                       -> do currenttandd <- DTime.getZonedTime
                                             _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                                    ++ "Removing downloaded BioMart region file ...")
                                             (_,_,_,ph) <- SP.createProcess $
                                                                (SP.shell ("rm "     ++
                                                                           outputdir ++
                                                                           "/"       ++
                                                                           "biomartresult.txt"
                                                                           ))
                                             ec <- waitForProcess ph
                                             case ec of
                                               SX.ExitFailure _ -> do currenttandd <- DTime.getZonedTime
                                                                      _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                                                             ++ "Could not remove downloaded BioMart region file ...")
                                                                      return processedregiondata
                                               SX.ExitSuccess   -> do currenttandd <- DTime.getZonedTime
                                                                      _ <- SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                                                             ++ "Removed downloaded BioMart region file ...")
                                                                      return processedregiondata
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
      outputdir     = DText.unpack $ outputdirectory config

{--------------------------}
