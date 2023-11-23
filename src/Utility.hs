{-# LANGUAGE MultiWayIf #-}

module Utility where

import Types

import qualified Data.ByteString                     as BS
import qualified Data.ByteString.UTF8                as UTF8
import           Data.CSV  as DCSV
import           Data.List as DL
import           Data.Text as DText
import           Data.Time as DTime
import           System.IO as SIO

numBytesUtf8Char :: Char
                 -> Int
numBytesUtf8Char =
  BS.length . UTF8.fromString . charToString
    where
      charToString :: Char -> String
      charToString c = [c]

numBytesUtf8String :: String
                   -> Int
numBytesUtf8String = BS.length . UTF8.fromString

stringToTuple :: [[String]]
              -> [[(String,String)]]
stringToTuple [] = []
stringToTuple xs = DL.map (DL.map (\y -> (y,"1")))
                          (DL.take ((DL.length xs) `div` 2) xs) 
                   DL.++ DL.map (DL.map (\y -> (y,"-1")))
                                (DL.drop ((DL.length xs) `div` 2) xs)

toCSV :: [[String]]
      -> String
toCSV allrows = genCsvFile allrows

writeCSV :: FRIConfig
         -> String
         -> String
         -> IO ()
writeCSV _      []      []             = return ()
writeCSV config csvdata outputfilename =
  if | DL.last outputdir == '/'
     -> SIO.writeFile (outputdir DL.++
                       outputfilename)
                       csvdata
     | otherwise
     -> SIO.writeFile (outputdir DL.++
                      "/"        DL.++ 
                      outputfilename)
                      csvdata
    where
      outputdir = DText.unpack $ outputdirectory config

showPrettyZonedTime :: ZonedTime
                    -> String
showPrettyZonedTime currenttandd =
  currenttanddnotimezone DL.++
  zeroestoadd            DL.++
  " "                    DL.++
  currenttimezone
    where
      zeroestoadd            = DL.concat                                               $
                               DL.map show                                             $
                               DL.take (26 - (DL.length currenttanddnotimezone))       $
                               DL.repeat 0
      currenttanddnotimezone = DL.reverse
                               (DL.drop 4
                               (DL.reverse
                               (show currenttandd)))
      currenttimezone        = DL.reverse
                               (DL.take 3
                               (DL.reverse
                               (show currenttandd)))
