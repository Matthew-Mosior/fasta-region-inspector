{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE KindSignatures    #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes        #-}
{-# LANGUAGE TypeOperators     #-}

module Linear.VariantLinear where

import AmbiguityCodes
import Variant
import Types

import           Data.Aeson.Types
import           Data.Text           as DText
import           Data.Vector.Unboxed as DVU
import           Effectful
import           Effectful.Ki
import           Effectful.Log

variantLinear :: forall {es :: [Effect]}.
                 ( StructuredConcurrency :> es
                 , Log :> es
                 , IOE :> es
                 )
              => [[(Variant,BioMartRegion,Char)]]
              -> [[[String]]]
              -> Eff es [[String]]
variantLinear xs
              ys =
  return $ variantsWithinAmbiguityCodesAndTSS xs
                                              ys
