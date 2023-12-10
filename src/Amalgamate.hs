{-# LANGUAGE MultiWayIf #-}

module Amalgamate where

import YamlParser

import Data.List as DL

amalgamateFinalVariantData :: [[String]]
                           -> [[String]]
                           -> [[String]]
amalgamateFinalVariantData []     [] = []
amalgamateFinalVariantData _      [] = []
amalgamateFinalVariantData []     _  = []
amalgamateFinalVariantData (x:xs) ys = 
  if | currentvariantidentifier `DL.elem` allfinalvariantidentifiers
     -> ((\(b:_) -> b) $
        DL.filter (\y -> (\(a:_) -> a) x == (\(a:_) -> a) y) ys)
        : amalgamateFinalVariantData xs ys
     | otherwise
     -> (DL.concat $
        ([x DL.++ ["N/A","N/A","N/A"]] :: [[String]]))
        : amalgamateFinalVariantData xs ys
    where
      currentvariantidentifier = (\(a:_) -> a)
                                 x
      allfinalvariantidentifiers = DL.map (\(a:_) -> a)
                                   ys
