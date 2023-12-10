{-# LANGUAGE MultiWayIf #-}

module Variant where

import Types
import YamlParser

import Data.List as DL
import Data.Text as DText

variantsAmbiguityCodesCheckerSmaller :: (Variant,BioMartRegion,Char)
                                     -> [String]
                                     -> Int
                                     -> [String]
variantsAmbiguityCodesCheckerSmaller _                                   []     _ = []
variantsAmbiguityCodesCheckerSmaller xs@(currentvariant,currentregion,_) (y:ys) z =
  --TSS reads in reverse direction (-1).
  if | currentregionstrand == "-1"
     -> if | (read y :: Int) >= (read currentvariantstartpos :: Int) &&
             (read currentvariantstartpos :: Int) >= ((((read y) - z) + 1) :: Int)
           -> [y] DL.++ (variantsAmbiguityCodesCheckerSmaller xs ys z)
           | otherwise
           -> variantsAmbiguityCodesCheckerSmaller xs ys z
     --TSS reads in forward direction (1).
     | (read y :: Int) <= (read currentvariantstartpos :: Int) &&
       (read currentvariantstartpos :: Int) <= ((((read y) + z) - 1) :: Int)
     -> [y] DL.++ (variantsAmbiguityCodesCheckerSmaller xs ys z)
     | otherwise
     -> variantsAmbiguityCodesCheckerSmaller xs ys z
    where
      currentregionstrand    = DText.unpack $ biomartregion_strand currentregion
      currentvariantstartpos = DText.unpack $ variant_startpos currentvariant      

variantsAmbiguityCodesCheckerSmall :: (Variant,BioMartRegion,Char)
                                   -> [[String]]
                                   -> [[String]]
variantsAmbiguityCodesCheckerSmall _                                               []     = []
variantsAmbiguityCodesCheckerSmall xs@(currentvariant,currentregion,withintssbool) (y:ys) =
  if | not (DL.null
           (variantsAmbiguityCodesCheckerSmaller xs
                                                 (DL.take ((DL.length y) - 6) (DL.drop 6 y))
                                                 (DL.length (y DL.!! 1))))
     -> [[DL.intercalate ":" allvariantdata]
         DL.++ [DL.intercalate ":" allregiondata]
         DL.++ [[withintssbool]]
         DL.++ [y DL.!! 0]
         DL.++ [y DL.!! 1]
         DL.++ [DL.intercalate ","
            (variantsAmbiguityCodesCheckerSmaller xs
            (DL.take ((DL.length y) - 6) (DL.drop 6 y))
            (DL.length (y DL.!! 1)))]]
         DL.++ (variantsAmbiguityCodesCheckerSmall xs ys)
     | otherwise
     -> variantsAmbiguityCodesCheckerSmall xs ys
    where
      allregiondata  = [ DText.unpack $ biomartregion_sequencedescription currentregion
                       , DText.unpack $ biomartregion_tss currentregion
                       , DText.unpack $ biomartregion_strand currentregion
                       , DText.unpack $ biomartregion_genename currentregion
                       ] 
      allvariantdata = [ DText.unpack $ variant_sample currentvariant
                       , DText.unpack $ variant_symbol currentvariant
                       , DText.unpack $ variant_sequencedescription currentvariant
                       , DText.unpack $ variant_startpos currentvariant
                       , DText.unpack $ variant_endpos currentvariant
                       , DText.unpack $ variant_ref currentvariant
                       , DText.unpack $ variant_alt currentvariant
                       , DText.unpack $ variant_enst currentvariant
                       ]

variantsAmbiguityCodesChecker :: (Variant,BioMartRegion,Char)
                              -> [[[String]]]
                              -> [[String]]
variantsAmbiguityCodesChecker _      []     = []
variantsAmbiguityCodesChecker xs     (y:ys) =
  (variantsAmbiguityCodesCheckerSmall xs y)
  DL.++ (variantsAmbiguityCodesChecker xs ys)

variantsWithinAmbiguityCodesAndTSSSmall :: [(Variant,BioMartRegion,Char)]
                                        -> [[[String]]]
                                        -> [[String]]
variantsWithinAmbiguityCodesAndTSSSmall []                          [] = []
variantsWithinAmbiguityCodesAndTSSSmall _                           [] = []
variantsWithinAmbiguityCodesAndTSSSmall []                          _  = []
variantsWithinAmbiguityCodesAndTSSSmall (x@(currentvariant,_,_):xs) ys =
  (variantsAmbiguityCodesChecker x ambiguitycodesregionsfiltered)
  DL.++ (variantsWithinAmbiguityCodesAndTSSSmall xs ys)
    where
      ambiguitycodesregionsfiltered = DL.map (DL.filter
                                             (\y -> (y DL.!! 5) ==
                                             currentvariantsymbol))
                                      ys
      currentvariantsymbol = DText.unpack $ variant_symbol currentvariant

variantsWithinAmbiguityCodesAndTSS :: [[(Variant,BioMartRegion,Char)]]
                                   -> [[[String]]]
                                   -> [[String]]
variantsWithinAmbiguityCodesAndTSS []     [] = []
variantsWithinAmbiguityCodesAndTSS _      [] = []
variantsWithinAmbiguityCodesAndTSS []     _  = []
variantsWithinAmbiguityCodesAndTSS (x:xs) ys =
  (variantsWithinAmbiguityCodesAndTSSSmall variantsfiltered ys) 
  DL.++ (variantsWithinAmbiguityCodesAndTSS xs ys)
    where 
      variantsfiltered = DL.filter (\(_,_,c) -> c == 'Y') x
