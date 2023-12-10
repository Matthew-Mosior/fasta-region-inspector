{-# LANGUAGE MultiWayIf #-}

module Region where

import Types

import Data.List  as DL
import Data.Maybe as DMaybe
import Data.Text  as DText

windowChecker :: Variant
              -> [BioMartRegion]
              -> Int
              -> [(Variant,BioMartRegion,Char)]
windowChecker _              []           _                           = []
windowChecker currentvariant (currentregion:restofregions) windowsize =
  if | currentregionstrand == "-1"      &&
       ((read currentregiontss :: Int) - windowsize) <=
       (read currentvariantstartpos :: Int) &&
       (read currentvariantstartpos :: Int) <=
       (read currentregiontss :: Int)
     -> [ ( currentvariant
          , currentregion
          , 'Y'
          )
        ]
        DL.++ (windowChecker currentvariant restofregions windowsize)
     | currentregionstrand == "1"       &&
       (read currentregiontss :: Int)       <=
       (read currentvariantstartpos :: Int) &&
       (read currentvariantstartpos :: Int) <=
       ((read currentregiontss :: Int) + windowsize)
     -> [ ( currentvariant
          , currentregion
          , 'Y'
          )
        ]
        DL.++ (windowChecker currentvariant restofregions windowsize)
     | otherwise
     -> [ ( currentvariant
          , currentregion
          , 'N'
          )
        ]
        DL.++ (windowChecker currentvariant restofregions windowsize)
    where
      currentregionstrand    = DText.unpack (biomartregion_strand currentregion)
      currentregiontss       = DText.unpack (biomartregion_tss currentregion)
      currentvariantstartpos = DText.unpack (variant_startpos currentvariant) 

variantRegionCheckSmall :: Variant
                        -> [BioMartRegion]
                        -> Int
                        -> [(Variant,BioMartRegion,Char)]
variantRegionCheckSmall _              []         _          = []
variantRegionCheckSmall currentvariant allregions windowsize =
  windowChecker currentvariant currentregion windowsize
  where
    currentregion = DL.filter (\x -> (biomartregion_genename x) ==
                                     (variant_symbol currentvariant)
                              )
                    allregions

variantRegionCheck :: [Variant]
                   -> [BioMartRegion]
                   -> Int
                   -> [[(Variant,BioMartRegion,Char)]]
variantRegionCheck _                               []         _          = []
variantRegionCheck []                              (_:_)      _          = []
variantRegionCheck (currentvariant:restofvariants) allregions windowsize =
  [variantRegionCheckSmall currentvariant allregions windowsize]
  DL.++ (variantRegionCheck restofvariants allregions windowsize)

variantWithinRegionCheck :: FRIConfig
                         -> [BioMartRegion]
                         -> [[(Variant,BioMartRegion,Char)]]
variantWithinRegionCheck _      []          = []
variantWithinRegionCheck config allregions  = 
  case DMaybe.isJust tsswinsize of
    True  -> do --User did provide tsswindowsize.
                variantRegionCheck (variants config)
                                   allregions
                                   (read $ DText.unpack $ DMaybe.fromJust tsswinsize)
    False -> do --User did not provide tsswindowsize.
                --Use 2 kb as default TSS window size.  
                variantRegionCheck (variants config)
                                   allregions
                                   2000
    where
      tsswinsize      = tsswindowsize config

convertToListWithinTSS :: [(Variant,BioMartRegion,Char)]
                       -> [String]
convertToListWithinTSS [] = []
convertToListWithinTSS xs = 
  DL.concat (DL.concat
            (DL.map (\(_,_,c) -> [[DL.intercalate ":" allvariantinfo] DL.++
                                 [DL.intercalate ":" allregioninfo]   DL.++
                                 [[c]]])
  xs))
    where
      allvariantinfo = DL.map (DText.unpack)
                              (DL.map (variant_sample)
                              (DL.map (\(x,_,_) -> x) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (variant_symbol)
                              (DL.map (\(x,_,_) -> x) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (variant_sequencedescription)
                              (DL.map (\(x,_,_) -> x) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (variant_startpos)
                              (DL.map (\(x,_,_) -> x) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (variant_endpos)
                              (DL.map (\(x,_,_) -> x) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (variant_ref)
                              (DL.map (\(x,_,_) -> x) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (variant_alt)
                              (DL.map (\(x,_,_) -> x) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (variant_enst)
                              (DL.map (\(x,_,_) -> x) xs))
      allregioninfo  = DL.map (DText.unpack)
                              (DL.map (biomartregion_sequencedescription)
                              (DL.map (\(_,y,_) -> y) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (biomartregion_tss)
                              (DL.map (\(_,y,_) -> y) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (biomartregion_strand)
                              (DL.map (\(_,y,_) -> y) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (biomartregion_genename)
                              (DL.map (\(_,y,_) -> y) xs))

prepareWithinTSS :: [[(Variant,BioMartRegion,Char)]]
                 -> [[String]]
prepareWithinTSS []     = []
prepareWithinTSS (x:xs) =
  [convertToListWithinTSS x]
  DL.++ (prepareWithinTSS xs)
