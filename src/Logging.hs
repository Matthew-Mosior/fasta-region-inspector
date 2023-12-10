{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE MultiWayIf    #-}
{-# LANGUAGE Strict        #-}

module Logging where

import Control.Concurrent (myThreadId,ThreadId)
import Data.List
import Data.Time
import GHC.Generics
import System.IO          (putStrLn)

data LoggingLevel = LogInfo
                  | LogDebug
                  | LogTrace
  deriving (Eq,Generic,Read,Show)

showPrettyLog :: LoggingLevel
              -> Int
              -> String
              -> String
              -> IO () 
showPrettyLog ll
              maxthreadnumber
              callingfunction
              message = do
  currenttandd                   <- getZonedTime 
  tid                            <- myThreadId
  let tidc                       = read   $ 
                                   drop 9 $ 
                                   show tid
  let currenttanddnotimezone     = reverse $
                                   drop 4  $
                                   reverse $
                                   show currenttandd
  let currenttimezone            = reverse $
                                   take 3  $
                                   reverse $
                                   show currenttandd
  let zeroestoadd                = concat                                                   $
                                   map show                                                 $
                                   take (zonedtimelength - (length currenttanddnotimezone)) $
                                   repeat 0
  let spacestoaddcallingfunction = take (longestfunctionnamelength - (length callingfunction)) $
                                   repeat ' '
  let spacestoaddtid             = if | tidc > 1 && tidc <= 9
                                      -> take 2 $
                                         repeat ' '
                                      | tidc >= 10 && tidc <= 99
                                      -> take 1 $
                                         repeat ' '
                                      | otherwise
                                      -> take 0 $
                                         repeat ' '
  putStrLn $ ( "["                        ++
               currenttanddnotimezone     ++
               zeroestoadd                ++
               " "                        ++
               currenttimezone            ++
               "]"                        ++
               " || "                     ++
               (show ll)                  ++
               " || "                     ++
               (show tid)                 ++
               spacestoaddtid             ++
               " || "                     ++
               callingfunction            ++
               spacestoaddcallingfunction ++
               " || "                     ++
               message
             )
    where
      zonedtimelength           = 29
      longestfunctionnamelength = 48
