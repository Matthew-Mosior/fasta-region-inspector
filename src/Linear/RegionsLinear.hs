{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE KindSignatures    #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes        #-}
{-# LANGUAGE TypeOperators     #-}

module Linear.RegionsLinear where

import AmbiguityCodes
import Region
import Types

import           Data.Aeson.Types
import           Data.Text           as DText
import           Data.Vector.Unboxed as DVU
import           Effectful
import           Effectful.Ki
import           Effectful.Log

regionsLinear :: forall {es :: [Effect]}.
                 ( StructuredConcurrency :> es
                 , Log :> es
                 , IOE :> es
                 )
              => Maybe Text
              -> FRIConfig
              -> [(String,String)]
              -> [[(String,String)]]
              -> [BioMartRegion]
              -> Eff es ([[(String,[String],[String],[[Int]])]],[[[String]]])
regionsLinear tsswinsizec
              config
              ambcodes
              mappedambiguitystrgroups
              allregions =
  case ignorestrandedness config of
    True -> do ambiguitycodeswithintss <- ambiguityCodesWithinRegionCheckIgnoreStrand tsswinsizec
                                                                                      config
                                                                                      ambcodes
                                                                                      mappedambiguitystrgroups
                                                                                      allregions
               --Prepare ambiguitycodeswithintss for printing.
               _ <- logMessage LogInfo
                               "Preparing ambiguity code strings final analysis." 
                               Null
               return $ (ambiguitycodeswithintss,prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss)
    False -> do ambiguitycodeswithintss <- ambiguityCodesWithinRegionCheck tsswinsizec
                                                                           config
                                                                           ambcodes
                                                                           mappedambiguitystrgroups
                                                                           allregions
                --Prepare ambiguitycodeswithintss for printing.
                _ <- logMessage LogInfo
                                "Preparing ambiguity code strings final analysis." 
                                Null 
                return $ (ambiguitycodeswithintss,prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss)
