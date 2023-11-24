{-# LANGUAGE DeriveGeneric #-}
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
              -> String
              -> String 
              -> IO () 
showPrettyLog ll callingfunction message = do
  currenttandd               <- getZonedTime 
  tid                        <- myThreadId
  let currenttanddnotimezone = reverse $
                               drop 4  $
                               reverse $
                               show currenttandd
  let currenttimezone        = reverse $
                               take 3  $
                               reverse $
                               show currenttandd
  let zeroestoadd            = concat                                      $
                               map show                                    $
                               take (26 - (length currenttanddnotimezone)) $
                               repeat 0
  putStrLn $ ( "["                    ++
               currenttanddnotimezone ++
               zeroestoadd            ++
               " "                    ++
               currenttimezone        ++
               "]"                    ++
               " || "                 ++
               (show ll)              ++
               " || "                 ++
               (show tid)             ++
               " || "                 ++
               callingfunction        ++
               " || "                 ++
               message
             )
