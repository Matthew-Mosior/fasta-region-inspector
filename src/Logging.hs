{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE DeriveGeneric     #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE RankNTypes        #-}
{-# LANGUAGE Strict            #-}
{-# LANGUAGE TypeFamilies      #-}
{-# LANGUAGE TypeOperators     #-}

module Logging where

import Control.Concurrent        (myThreadId,ThreadId)
import Data.List
import Data.Time
import Effectful
import Effectful.Dispatch.Static
import GHC.Generics
import System.IO                 (putStrLn)

data LoggingLevel = LogInfo
                  | LogDebug
                  | LogTrace
  deriving (Eq,Generic,Read,Show)

data FRILogging :: Effect

type instance DispatchOf FRILogging = 'Static 'WithSideEffects
data instance StaticRep  FRILogging = FRILogging 

runFRILogging :: IOE :> es => Eff (FRILogging : es) a -> Eff es a
runFRILogging = evalStaticRep FRILogging

showPrettyLog :: forall {es :: [Effect]} {b}.
                 ( FRILogging :> es
                 , IOE :> es
                 )
              => LoggingLevel
              -> Int
              -> String
              -> String
              -> Eff es () 
showPrettyLog ll
              maxthreadnumber
              callingfunction
              message = do
  currenttandd                   <- liftIO getZonedTime 
  tid                            <- liftIO myThreadId
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
  liftIO     $
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
