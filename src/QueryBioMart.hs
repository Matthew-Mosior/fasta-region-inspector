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

import Control.Monad as CM
import Control.Monad.IO.Class as CMIOC
import Control.Exception as CE
import Control.Retry as CR
import Data.ByteString.Char8 as DBC8
import Data.Function as DF
import Data.List as DL
import Data.List.Split as DLS
import Data.Text as DText
import Data.Text.Lazy as DTextL
import qualified Data.Text.Lazy as DTL
import qualified Data.Text.Lazy.IO as DTLIO
import qualified Data.Map as Map
import Data.Time as DTime
import Network.Connection as NC
import Network.HTTP.Client as NHTTPC
import Network.HTTP.Client.TLS as NHTTPTLS
import Network.HTTP.Types as NHTTPT
import Network.HTTP.Req as NHR
import Network.TLS as NTLS
import Network.TLS.Extra.Cipher as NTLSEC
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
  --Create new http manager.
  mymanager <- NHTTPC.newManager $ mkManagerSettings (TLSSettingsSimple True False False) Nothing
  --Define httpconfig using mymanager.
  let httpconfig = HttpConfig { httpConfigProxy = Nothing
                              , httpConfigRedirectCount = 10
                              , httpConfigAltManager = Just mymanager
                              , httpConfigCheckResponse = \_ response preview ->
                                  let scode = statusCodeR response
                                    in if | 200 <= scode && scode < 300
                                          -> Nothing
                                          | otherwise
                                          -> Just (StatusCodeException (void response) preview)
                              , httpConfigRetryPolicy = retryPolicyDefault
                              , httpConfigRetryJudge  = \_ response ->
                                  statusCodeR response
                                    `DL.elem` ([ 408 -- Request timeout
                                              , 504 -- Gateway timeout
                                              , 524 -- A timeout occured
                                              , 598 -- (Informal convention) Network read timeout error
                                              , 599 -- (Informal convention) Network connect timeout error
                                              ] :: [Int])
                              , httpConfigRetryJudgeException = \_ e ->
                                  case fromException e of
                                    Just (HttpExceptionRequest _ c) ->
                                      case c of
                                        ResponseTimeout -> True
                                        ConnectionTimeout -> True
                                        _ -> False
                                    _ -> False
                              , httpConfigBodyPreviewLength = 1024
                              }
                                where
                                  statusCodeR = NHTTPT.statusCode DF.. responseStatus
  runReq httpconfig $ do
    currenttandd <- liftIO DTime.getZonedTime
    _ <- liftIO $ SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                    ++ "Querying and downloading region data from BioMart via HTTP request ...")
    let params = "query" =: ((DTextL.toStrict $ renderText def finalbiomartxml) :: DText.Text)
    biomartrequest <- req NHR.POST (http "www.ensembl.org" /: "biomart" /: "martservice")
                                   (ReqBodyUrlEnc params)
                                   bsResponse
                                   mempty
    currenttandd <- liftIO DTime.getZonedTime
    _ <- liftIO $ SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                    ++ "Successfully queried and returned region data from BioMart via HTTP request ...")
    let returnedbiomartregions = DBC8.unpack $
                                 NHR.responseBody biomartrequest 
    let finalbiomartregions = DL.map (\x -> DLS.splitOn "," x)
                              (DL.lines returnedbiomartregions)
    let processedregiondata = DL.map (\x -> BioMartRegion { rchromosome = DText.pack
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
       -> if | DL.last outputdir == '/'
             -> do currenttandd <- liftIO DTime.getZonedTime
                   _ <- liftIO $ SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                   ++ "Writing BioMart region data to file biomartresult.txt in output directory ...")
                   _ <- liftIO $ SIO.writeFile (outputdir ++ "biomartresult.txt")
                                               returnedbiomartregions
                   return processedregiondata
             | otherwise
             -> do currenttandd <- liftIO DTime.getZonedTime
                   _ <- liftIO $ SIO.putStrLn ("[" ++ (showPrettyZonedTime currenttandd) ++ "] "
                                                   ++ "Writing BioMart region data to file biomartresult.txt in output directory ...")
                   _ <- liftIO $ SIO.writeFile (outputdir ++ "/" ++ "biomartresult.txt")
                                               returnedbiomartregions
                   return processedregiondata 
       | otherwise
       -> return processedregiondata
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
