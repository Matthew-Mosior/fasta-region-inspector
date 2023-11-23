{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE KindSignatures    #-}
{-# LANGUAGE MultiWayIf        #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes        #-}
{-# LANGUAGE TypeOperators     #-}

module AmbiguityCodes where

import Types

import Linear.UtilityLinear

import Codec.Binary.UTF8.String as CBUTF8
import Control.Parallel.Strategies as CPS
import Data.Aeson.Types
import Data.ByteString as DB hiding (append)
import Data.ByteString.Char8 as DBC hiding (append)
import Data.ByteString.Search.DFA as DBSDFA
import Data.Char as DC
import Data.Either as DE
import Data.List as DL
import Data.List.Split as DLS
import Data.Maybe as DMaybe
import Data.SBV as DSBV
import qualified Data.SBV.String as DSBVS
import Data.SBV.RegExp as DSBVRE
import Data.Text as DText
import Data.Traversable as DT
import Data.Vector.Storable.ByteString as DVSBS
import Data.Vector.Unboxed as DVU
import Effectful
import Effectful.Ki
import Effectful.Log

reverseComplementNucleotide :: DB.ByteString
                            -> DB.ByteString
reverseComplementNucleotide currentsequence =
  DB.pack
  (CBUTF8.encode (DL.map snd
                 (DL.concatMap
                 (\y -> DL.filter (\(r,_) -> r == y) revcomplementmapping)
                 (CBUTF8.decode
                 (DB.unpack (DB.reverse currentsequence))))))
    where
      revcomplementmapping = [ ('A','T')
                             , ('T','A')
                             , ('G','C')
                             , ('C','G')
                             ]

grabRegionSequence :: Vector Char
                   -> BioMartRegion
                   -> Int
                   -> DB.ByteString
grabRegionSequence currentsequence
                   currentregion
                   currenttsswinsize =
  if | (currentregionstrand == "-1")
     -> DB.pack
        (CBUTF8.encode
        (DVU.toList
        (DVU.take currenttsswinsize (DVU.drop ((read currentregiontss) - (currenttsswinsize - 1))
                                    currentsequence))))
     | otherwise
     -> DB.pack
        (CBUTF8.encode
        (DVU.toList
        (DVU.take currenttsswinsize (DVU.drop ((read currentregiontss) - 1)
                                   currentsequence))))
    where
      currentregiontss    = DText.unpack (biomartregion_tss currentregion)
      currentregionstrand = DText.unpack (biomartregion_strand currentregion)

subStrLocationsSmallForward :: ( MonadIO m
                               , MonadLog m
                               )
                            => Maybe Text
                            -> [String]
                            -> BioMartRegion
                            -> Vector Char
                            -> m [[Int]]
subStrLocationsSmallForward _ [] _ _ = return []
subStrLocationsSmallForward tsswinsizec
                            (currentmappedambstr:restofmappedambstrs)
                            currentregion
                            finalfastafile = do
  _ <- logMessage LogInfo
                  ( "Processing mapped ambiguity code "
                    `append`
                    (DText.pack currentmappedambstr)
                    `append`
                    "."
                  )
                  Null
  if | DMaybe.isJust tsswinsizec
     -> do res <- subStrLocationsSmallForward tsswinsizec
                                              restofmappedambstrs
                                              currentregion
                                              finalfastafile
           return $ (DBSDFA.indices (DBC.pack currentmappedambstr)
                                    (grabRegionSequence
                                    finalfastafile
                                    currentregion
                                    (read (DText.unpack $ DMaybe.fromJust tsswinsizec) :: Int))) : res
     | otherwise
     -> do res <- subStrLocationsSmallForward tsswinsizec
                                              restofmappedambstrs
                                              currentregion
                                              finalfastafile
           return $ (DBSDFA.indices (DBC.pack currentmappedambstr)
                                    (grabRegionSequence
                                    finalfastafile
                                    currentregion
                                    2000)) : res

subStrLocationsSmallReverse :: ( MonadIO m
                               , MonadLog m
                               )
                            => Maybe Text
                            -> [String]
                            -> BioMartRegion
                            -> Vector Char
                            -> m [[Int]]
subStrLocationsSmallReverse _ [] _ _ = return []
subStrLocationsSmallReverse tsswinsizec
                            (currentmappedambstr:restofmappedambstrs)
                            currentregion
                            finalfastafile = do
  _ <- logMessage LogInfo
                  ( "Processing mapped ambiguity code "
                    `append`
                    (DText.pack currentmappedambstr)
                    `append`
                    "."
                  )
                  Null
  if | DMaybe.isJust tsswinsizec
     -> do res <- subStrLocationsSmallReverse tsswinsizec
                                              restofmappedambstrs
                                              currentregion
                                              finalfastafile
           return $ (DL.map (\a -> (DBC.length
                                   ((grabRegionSequence finalfastafile
                                                        currentregion
                                                        (read (DText.unpack $ DMaybe.fromJust tsswinsizec) :: Int)))) - a - 1)
                                   (DBSDFA.indices (DBC.pack currentmappedambstr)
                                   (reverseComplementNucleotide
                                   (grabRegionSequence finalfastafile
                                                       currentregion
                                                       (read (DText.unpack $ DMaybe.fromJust tsswinsizec) :: Int))))) : res
     | otherwise
     -> do res <- subStrLocationsSmallReverse tsswinsizec
                                              restofmappedambstrs
                                              currentregion
                                              finalfastafile
           return $ (DL.map (\a -> (DBC.length
                                   ((grabRegionSequence finalfastafile
                                                        currentregion
                                                        2000))) - a - 1)
                                   (DBSDFA.indices (DBC.pack currentmappedambstr)
                                   (reverseComplementNucleotide
                                   (grabRegionSequence finalfastafile
                                                       currentregion
                                                       2000)))) : res

subStrLocations :: ( MonadIO m
                   , MonadLog m
                   )
                => Maybe Text
                -> [String]
                -> BioMartRegion
                -> Vector Char
                -> m [[Int]]
subStrLocations _ [] _ _ = return []
subStrLocations tsswinsizec
                allmappedambiguitystrs
                currentregion
                fastaseq = do
  if | DMaybe.isJust tsswinsizec
     -> if | currentregionstrand == "-1"
           -> do reversesubstrlocs <- subStrLocationsSmallReverse tsswinsizec
                                                                  allmappedambiguitystrs
                                                                  currentregion
                                                                  fastaseq
                 return $ ((DL.map (DL.map (\i ->
                            ((((read currentregiontss) - (read (DText.unpack $ DMaybe.fromJust tsswinsizec) :: Int)) + i) + 2)))
                          reversesubstrlocs)
                          `CPS.using` (CPS.parList CPS.rdeepseq))
           | otherwise
           -> do forwardsubstrlocs <- subStrLocationsSmallForward tsswinsizec
                                                                  allmappedambiguitystrs
                                                                  currentregion
                                                                  fastaseq
                 return $ ((DL.map (DL.map (\i ->
                            ((read currentregiontss) + i)))
                          forwardsubstrlocs)
                          `CPS.using` (CPS.parList CPS.rdeepseq))
     | otherwise
     -> if | currentregionstrand == "-1"
           -> do reversesubstrlocs <- subStrLocationsSmallReverse tsswinsizec
                                                                  allmappedambiguitystrs
                                                                  currentregion
                                                                  fastaseq
                 return $ ((DL.map (DL.map (\i ->
                            ((((read currentregiontss) - 2000) + i) + 2)))
                          reversesubstrlocs)
                          `CPS.using` (CPS.parList CPS.rdeepseq))
           | otherwise
           -> do forwardsubstrlocs <- subStrLocationsSmallForward tsswinsizec
                                                                  allmappedambiguitystrs
                                                                  currentregion
                                                                  fastaseq
                 return $ ((DL.map (DL.map (\i ->
                            ((read currentregiontss) + i)))
                          forwardsubstrlocs)
                          `CPS.using` (CPS.parList CPS.rdeepseq))
    where
      currentregionstrand = DText.unpack $ biomartregion_strand currentregion
      currentregiontss    = DText.unpack $ biomartregion_tss currentregion

ambiguityCodesWithinRegionCheckIgnoreStrandSmall :: forall {es :: [Effect]} {b}.
                                                    ( StructuredConcurrency :> es
                                                    , Log :> es
                                                    , IOE :> es
                                                    )
                                                 => Maybe Text
                                                 -> FRIConfig
                                                 -> (String,String)
                                                 -> [(String,String)]
                                                 -> [BioMartRegion]
                                                 -> Eff es [(String,[String],[String],[[Int]])]
ambiguityCodesWithinRegionCheckIgnoreStrandSmall _ _ ([],[])       []    _  = return []
ambiguityCodesWithinRegionCheckIgnoreStrandSmall _ _ _             []    _  = return []
ambiguityCodesWithinRegionCheckIgnoreStrandSmall _ _ ([],[])       _     _  = return []
ambiguityCodesWithinRegionCheckIgnoreStrandSmall _ _ ((_:_),[])    (_:_) [] = return []
ambiguityCodesWithinRegionCheckIgnoreStrandSmall _ _ ((_:_),(_:_)) (_:_) [] = return []
ambiguityCodesWithinRegionCheckIgnoreStrandSmall _ _ ([],(_:_))    (_:_) [] = return []
ambiguityCodesWithinRegionCheckIgnoreStrandSmall tsswinsizec
                                                 config
                                                 currentambtuple@(currentambcode,currentstrand)
                                                 allmappedambiguitystrs
                                                 allregions =
  scoped $ \scope -> do
    DT.forM allregions $ \currentregion -> do
      regioncheckdata <- fork scope (do let currentregionchr      = DText.unpack (biomartregion_sequencedescription currentregion)
                                        let currentregiontss      = DText.unpack (biomartregion_tss currentregion)
                                        let currentregionstrand   = DText.unpack (biomartregion_strand currentregion)
                                        let currentregiongenename = DText.unpack (biomartregion_genename currentregion)
                                        let allcurrentregiondata  = [ currentregionchr
                                                                    , currentregiontss
                                                                    , currentregionstrand
                                                                    , currentregiongenename
                                                                    ]
                                        _ <- logMessage LogInfo
                                                        ( "Processing region data associated with gene "
                                                          `append`
                                                          (DText.pack currentregiongenename)
                                                          `append`
                                                          "."
                                                        )
                                                        Null
                                        --Ignore strandedness, find both the ambiguity mapped strings
                                        --and the reverse complement ambiguity mapped strings.
                                        --Grab locations of mapped ambiguity codes,
                                        --and recurse.
                                        --Grab locations of mapped am codes,
                                        --and recurse.
                                        cfai       <- liftIO $ getFAILineLinear config
                                                                                currentregion
                                                                                (toInteger 0)
                                        fastaseq   <- liftIO $ getFASTASequenceLinear config
                                                                                      cfai
                                                                                      (toInteger 0)
                                                                                      DVU.empty 
                                        substrlocs <- subStrLocations tsswinsizec
                                                                      (DL.map fst allmappedambiguitystrs)
                                                                      currentregion
                                                                      ((\(FASTASequence seq) -> seq) fastaseq)
                                        return $ ( currentambcode
                                                 , DL.map fst allmappedambiguitystrs
                                                 , allcurrentregiondata
                                                 , substrlocs
                                                 )
                                    )
      atomically $ await regioncheckdata

ambiguityCodesWithinRegionCheckSmall :: forall {es :: [Effect]} {b}.
                                        ( StructuredConcurrency :> es
                                        , Log :> es
                                        , IOE :> es
                                        )
                                     => Maybe Text
                                     -> FRIConfig
                                     -> ([Char], [Char])
                                     -> [(String, b)]
                                     -> [BioMartRegion]
                                     -> Eff es [([Char], [String], [String], [[Int]])]
ambiguityCodesWithinRegionCheckSmall _ _ ([],[])       []    _  = return []
ambiguityCodesWithinRegionCheckSmall _ _ _             []    _  = return []
ambiguityCodesWithinRegionCheckSmall _ _ ([],[])       _     _  = return []
ambiguityCodesWithinRegionCheckSmall _ _ ((_:_),[])    (_:_) [] = return []
ambiguityCodesWithinRegionCheckSmall _ _ ((_:_),(_:_)) (_:_) [] = return []
ambiguityCodesWithinRegionCheckSmall _ _ ([],(_:_))    (_:_) [] = return []
ambiguityCodesWithinRegionCheckSmall tsswinsizec
                                     config
                                     currentambtuple@(currentambcode,currentstrand)
                                     allmappedambiguitystrs
                                     allregions =
  scoped $ \scope -> do
    DT.forM allregions $ \currentregion -> do
      regioncheckdata <- fork scope (do let currentregionchr      = DText.unpack (biomartregion_sequencedescription currentregion)
                                        let currentregiontss      = DText.unpack (biomartregion_tss currentregion)
                                        let currentregionstrand   = DText.unpack (biomartregion_strand currentregion)
                                        let currentregiongenename = DText.unpack (biomartregion_genename currentregion)
                                        let allcurrentregiondata  = [ currentregionchr
                                                                    , currentregiontss
                                                                    , currentregionstrand
                                                                    , currentregiongenename
                                                                    ]
                                        _ <- logMessage LogInfo
                                                        ( "Processing region data associated with gene "
                                                          `append`
                                                          (DText.pack currentregiongenename)
                                                          `append`
                                                          "."
                                                        )
                                                        Null
                                        if | currentregionstrand == currentstrand
                                           -> do --Grab locations of mapped am codes,
                                                 --and recurse.
                                                 cfai       <- liftIO $ getFAILineLinear config
                                                                                         currentregion
                                                                                         (toInteger 0)
                                                 fastaseq   <- liftIO $ getFASTASequenceLinear config
                                                                                               cfai
                                                                                               (toInteger 0)
                                                                                               DVU.empty 
                                                 substrlocs <- subStrLocations tsswinsizec
                                                                               (DL.map fst allmappedambiguitystrs)
                                                                               currentregion
                                                                               ((\(FASTASequence seq) -> seq) fastaseq)
                                                 return $ ( currentambcode
                                                          , DL.map fst allmappedambiguitystrs
                                                          , allcurrentregiondata
                                                          , substrlocs
                                                          )
                                           | otherwise
                                           -> do --Current ambiguity codes and mapped strings
                                                 --are not correct for region strand
                                                 --(i.e. "-1" != "1" or "1" != "-1").
                                                 let numofspaces      = DText.length "                               "
                                                 let printnumofspaces = DText.replicate numofspaces (DText.singleton ' ')
                                                 _ <- logMessage LogInfo
                                                                 ( "Could not process region data associated with current ambiguity code "
                                                                   `append`
                                                                   (DText.pack currentambcode)
                                                                   `append`
                                                                   ":\n"
                                                                   `append`
                                                                   printnumofspaces
                                                                   `append`
                                                                   (DText.pack currentambcode)
                                                                   `append`
                                                                   " strand orientation is "
                                                                   `append`
                                                                   (DText.pack currentregionstrand)
                                                                   `append`
                                                                   "."
                                                                 )
                                                                 Null
                                                 return $ ( currentambcode
                                                          , DL.map fst allmappedambiguitystrs
                                                          , allcurrentregiondata
                                                          , []
                                                          )
                                    )
      atomically $ await regioncheckdata                     

ambiguityCodesWithinRegionCheckIgnoreStrand :: forall {es :: [Effect]} {b}.
                                               ( StructuredConcurrency :> es
                                               , Log :> es
                                               , IOE :> es
                                               )
                                            => Maybe Text
                                            -> FRIConfig
                                            -> [(String,String)]
                                            -> [[(String,String)]]
                                            -> [BioMartRegion]
                                            -> Eff es [[(String,[String],[String],[[Int]])]]
ambiguityCodesWithinRegionCheckIgnoreStrand _ _ [] [] _ = return []
ambiguityCodesWithinRegionCheckIgnoreStrand _ _ _  [] _ = return []
ambiguityCodesWithinRegionCheckIgnoreStrand _ _ [] _  _ = return []
ambiguityCodesWithinRegionCheckIgnoreStrand tsswinsizec
                                            config
                                            (currentambcode:restofambcodes)
                                            (currentmappedambiguitystrgroup:restofmappedambiguitystrgroups)
                                            allregions = do
  ambcodeswithinregioncheck <- ambiguityCodesWithinRegionCheckIgnoreStrandSmall tsswinsizec
                                                                                config
                                                                                currentambcode
                                                                                currentmappedambiguitystrgroup
                                                                                allregions
  res <- ambiguityCodesWithinRegionCheckIgnoreStrand tsswinsizec
                                                     config
                                                     restofambcodes
                                                     restofmappedambiguitystrgroups
                                                     allregions
  return $ ambcodeswithinregioncheck : res

ambiguityCodesWithinRegionCheck :: forall {es :: [Effect]} {b}.
                                   ( StructuredConcurrency :> es
                                   , Log :> es
                                   , IOE :> es
                                   )
                                => Maybe Text
                                -> FRIConfig
                                -> [(String,String)]
                                -> [[(String,String)]]
                                -> [BioMartRegion]
                                -> Eff es [[(String,[String],[String],[[Int]])]]
ambiguityCodesWithinRegionCheck _ _ [] [] _ = return []
ambiguityCodesWithinRegionCheck _ _ _  [] _ = return []
ambiguityCodesWithinRegionCheck _ _ [] _  _ = return []
ambiguityCodesWithinRegionCheck tsswinsizec
                                config
                                (currentambcode:restofambcodes)
                                (currentmappedambiguitystrgroup:restofmappedambiguitystrgroups)
                                allregions = do
  ambcodeswithinregioncheck <- ambiguityCodesWithinRegionCheckSmall tsswinsizec
                                                                    config
                                                                    currentambcode
                                                                    currentmappedambiguitystrgroup
                                                                    allregions
  res <- ambiguityCodesWithinRegionCheck tsswinsizec
                                         config
                                         restofambcodes
                                         restofmappedambiguitystrgroups
                                         allregions
  return $ ambcodeswithinregioncheck : res
 
allStrGenerationSmall :: String
                      -> IO [String]
allStrGenerationSmall [] = return []
allStrGenerationSmall xs = do 
  --Need to create all possible strings created from nucleotidemapping.
  --Use Data.SBV library to create all possible satisfiability predicates
  --(Strings generated from a particular regular expression pattern).
  generatedregexstrs <- DSBV.allSat $ \s -> (s :: SString) 
                                            `DSBVRE.match` 
                                            (DSBVRE.Conc 
                                            (DL.map (DSBVRE.oneOf)
                                            nucleotidemappingfinal))
  --Return only true results of generatedregexstrs.        
  return (DL.filter (\z -> not (DL.null z)) 
                    (DL.filter (\y -> DL.all (DC.isUpper) y) 
                    (DLS.splitOneOf " " 
                    (DL.filter (\x -> not (x `DL.elem` ("\"" :: String))) 
                    (show generatedregexstrs)))))
    where
      charConversion :: String -> [(Char,Char)] -> [(Char,String)]
      charConversion []     [] = []
      charConversion _      [] = []
      charConversion []     _  = []
      charConversion (x:xs) ys = [(DL.head (DL.map (fst) 
                                                   (DL.filter (\(a,_) -> a == x) ys)),
                                                   DL.map (snd) (DL.filter (\(a,_) -> a == x) ys))] 
                                 DL.++ (charConversion xs ys)
      nucleotidemappingfinal    = DL.map (snd) 
                                  (charConversion xs nucleotidemappingfiltered) 
      nucleotidemappingfiltered = DL.filter (\(b,_) -> b 
                                                       `DL.elem` 
                                                       (DL.filter 
                                                       (\a -> a `DL.elem` DL.map (fst) nucleotidemapping) xs)) 
                                                       nucleotidemapping 
      nucleotidemapping = [ ('A','A')
                          , ('T','T')
                          , ('G','G')
                          , ('C','C')
                          , ('Y','C')
                          , ('Y','T')
                          , ('R','A')
                          , ('R','G')
                          , ('W','A')
                          , ('W','T')
                          , ('S','G')
                          , ('S','C')
                          , ('K','T')
                          , ('K','G')
                          , ('M','C')
                          , ('M','A')
                          , ('D','A')
                          , ('D','G')
                          , ('D','T')
                          , ('V','A')
                          , ('V','C')
                          , ('V','G')
                          , ('H','A')
                          , ('H','C')
                          , ('H','T')
                          , ('B','C')
                          , ('B','G')
                          , ('B','T')
                          , ('N','A')
                          , ('N','T')
                          , ('N','G')
                          , ('N','C')
                          , ('X','A')
                          , ('X','T')
                          , ('X','G')
                          , ('X','C')
                          ]

allStrGeneration :: [String]
                 -> IO [[String]]
allStrGeneration [] = return []
allStrGeneration xs = DT.mapM (allStrGenerationSmall) xs

ambiguityCodesReverseComplementSmall :: [String]
                                     -> [String]
ambiguityCodesReverseComplementSmall []                               = []
ambiguityCodesReverseComplementSmall (currentambcodes:restofambcodes) = 
  [reversecomplementfinal]
  DL.++ (ambiguityCodesReverseComplementSmall restofambcodes)
    where
      reversecomplementfinal = DL.map (snd)
                                      (DL.concatMap (\y -> DL.filter (\(r,_) -> r == y)
                                      reversecomplementfiltered)
                                      (DL.reverse currentambcodes))
      reversecomplementfiltered = DL.filter (\(b,_) -> b
                                                       `DL.elem`
                                                       (DL.filter (\a -> a
                                                                         `DL.elem`
                                                                         DL.map
                                                                         (fst) revcomplementmapping)
                                                       currentambcodes))
                                            revcomplementmapping
      revcomplementmapping = [ ('A','T')
                             , ('T','A')
                             , ('G','C')
                             , ('C','G')
                             , ('Y','R')
                             , ('R','Y')
                             , ('W','W')
                             , ('S','S')
                             , ('K','M')
                             , ('M','K')
                             , ('D','H')
                             , ('H','D')
                             , ('V','B')
                             , ('B','V')
                             , ('X','X')
                             , ('N','N')
                             , ('-','-')
                             ]

ambiguityCodesReverseComplement :: FRIConfig
                                -> [String]
ambiguityCodesReverseComplement config =
  ambiguityCodesReverseComplementSmall ambcodes  
    where
      ambcodes = DL.map DText.unpack
                 (ambiguitycodes config)

tupleConverterAmbCodesWithinTSS :: (String,[String],[String],[[Int]])
                                -> ([[String]],[[String]],[[String]],[[Int]])
tupleConverterAmbCodesWithinTSS (a,b,c,d) =
  ((DL.map (\x -> [x])
   (DL.replicate (DL.length b) a))
  ,(DL.map (\x -> [x]) b)
  ,(DL.replicate (DL.length b) c)
  ,d)

tupleToListAmbCodesWithinTSS :: ([[String]],[[String]],[[String]],[[Int]])
                             -> [[String]]
tupleToListAmbCodesWithinTSS ([],    _,     _,     _)      = []
tupleToListAmbCodesWithinTSS ((_:_), [],    _,     _)      = []
tupleToListAmbCodesWithinTSS ((_:_), (_:_), [],    _)      = []
tupleToListAmbCodesWithinTSS ((_:_), (_:_), (_:_), [])     = []
tupleToListAmbCodesWithinTSS ((a:as),(b:bs),(c:cs),(d:ds)) =
  [a DL.++ b DL.++ c DL.++ (DL.map (show) d)]
  DL.++ (tupleToListAmbCodesWithinTSS (as,bs,cs,ds))

convertToListAmbCodesWithinTSS :: [(String,[String],[String],[[Int]])]
                               -> [[[String]]]
convertToListAmbCodesWithinTSS [] = []
convertToListAmbCodesWithinTSS xs =
  DL.map (\y -> tupleToListAmbCodesWithinTSS y)
         (DL.map (\x -> tupleConverterAmbCodesWithinTSS x)
  xs)

prepareAmbiguityCodesWithinTSS :: [[(String,[String],[String],[[Int]])]]
                               -> [[[String]]]
prepareAmbiguityCodesWithinTSS []     = []
prepareAmbiguityCodesWithinTSS (x:xs) =
  (convertToListAmbCodesWithinTSS x)
  DL.++ (prepareAmbiguityCodesWithinTSS xs)
