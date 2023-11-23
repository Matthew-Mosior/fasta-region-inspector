{-# LANGUAGE OverloadedStrings #-}

module Linear.VariantLinear where

import AmbiguityCodes
import Variant
import Types

import           Data.Aeson.Types
import           Data.Text           as DText
import           Data.Vector.Unboxed as DVU
import           Effectful
import           Effectful.Log

variantLinear :: ( MonadIO m
                 , MonadLog m
                 )
              => [[(Variant,BioMartRegion,Char)]]
              -> [[[String]]]
              -> m [[String]]
variantLinear xs
              ys =
  return $ variantsWithinAmbiguityCodesAndTSS xs
                                              ys
