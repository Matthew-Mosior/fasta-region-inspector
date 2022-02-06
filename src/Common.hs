{-=Fasta-Region-Inspector (FRI): A Somatic=-}
{-=Hypermutation Analysis Tool.=-}
{-=Author: Matthew Mosior=-}
{-=Synposis: This Haskell script will=-}
{-=process command line arguments to FRI.=-}


{-Language extension.-}

{-# LANGUAGE Strict     #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE MultiWayIf #-}

{---------------------}


{-Module.-}

module Common where

{---------}


{-Import modules.-}

import YamlParser

{-----------------}


{-Imports.-}

import Codec.Binary.UTF8.String as CBUTF8
import Control.DeepSeq as CD
import Control.Parallel.Strategies as CPS
import Data.Attoparsec.ByteString as DAB
import Data.ByteString as DB
import Data.ByteString.Char8 as DBC
import Data.ByteString.Lazy as DBL
import Data.ByteString.Search.DFA as DBSDFA
import Data.Char as DC
import Data.CSV as DCSV
import Data.Either as DE
import Data.Functor as DF
import Data.List as DL
import Data.List.Split as DLS
import Data.Maybe as DMaybe
import Data.Ord as DO
import Data.SBV as DSBV
import qualified Data.SBV.String as DSBVS
import Data.SBV.RegExp as DSBVRE
import Data.Text as DText
import Data.Time as DTime
import Data.Traversable as DT
import Data.Vector.Storable.ByteString as DVSBS
import Data.Vector.Unboxed as DVU
import ELynx.Import.Sequence.Fasta as EISF
import qualified ELynx.Data.Alphabet.Alphabet as EDAA
import qualified ELynx.Data.Alphabet.Character as EDAC
import qualified ELynx.Data.Character.Character as EDCC
import qualified ELynx.Data.Character.NucleotideI as EDCN
import qualified ELynx.Data.Sequence.Sequence as EDSS
import System.Console.GetOpt as SCG
import System.Process as SP
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO
import System.IO.Temp as SIOT
import Text.Regex.TDFA as TRTDFA

{----------}


{-General Utility Functions.-}

--stringToTuple -> This function will
--annotate all mapped strings with
--directionality (tuple). 
stringToTuple :: [[String]] -> [[(String,String)]]
stringToTuple [] = []
stringToTuple xs = DL.map (DL.map (\y -> (y,"1")))
                          (DL.take ((DL.length xs) `div` 2) xs) 
                   DL.++ DL.map (DL.map (\y -> (y,"-1")))
                                (DL.drop ((DL.length xs) `div` 2) xs)

{----------------------------}

{-Region functions.-}

--windowChecker -> This function will
--check to see if the given variant is within genes
--2 kb (or custom TSS window size if provided)
--of the genes TSS.
windowChecker :: Variants -> [BioMartRegion] -> Int -> [(Variants,BioMartRegion,Char)]
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
      currentregionstrand    = DText.unpack (rstrand currentregion)
      currentregiontss       = DText.unpack (rtss currentregion)
      currentvariantstartpos = DText.unpack (vstartpos currentvariant) 

--variantRegionCheckSmall -> This function will
--check to see if the given variant is within genes
--2 kb (or custom TSS window size if provided)
--of the genes TSS.
variantRegionCheckSmall :: Variants -> [BioMartRegion] -> Int -> [(Variants,BioMartRegion,Char)]
variantRegionCheckSmall _              []         _          = []
variantRegionCheckSmall currentvariant allregions windowsize =
  windowChecker currentvariant currentregion windowsize
  where
    currentregion = DL.filter (\x -> (rgenename x) ==
                                     (vsymbol currentvariant)
                              )
                    allregions

--variantRegionCheck -> This function will
--check to see if the given variant is within genes
--2 kb (or custom TSS window size if provided)
--of the genes TSS.
variantRegionCheck :: [Variants] -> [BioMartRegion] -> Int -> [[(Variants,BioMartRegion,Char)]]
variantRegionCheck _                               []         _          = []
variantRegionCheck []                              (_:_)      _          = []
variantRegionCheck (currentvariant:restofvariants) allregions windowsize =
  [variantRegionCheckSmall currentvariant allregions windowsize]
  DL.++ (variantRegionCheck restofvariants allregions windowsize)

--variantWithinRegionCheck -> This function will
--check to see if the given variant is within genes
--2 kb (or custom TSS window size if provided)
--of the genes TSS.
variantWithinRegionCheck :: FRIConfig -> [BioMartRegion] -> [[(Variants,BioMartRegion,Char)]]
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

--convertToListWithinTSS -> This function will
--prepare for printFile function.
convertToListWithinTSS :: [(Variants,BioMartRegion,Char)] -> [String]
convertToListWithinTSS [] = []
convertToListWithinTSS xs = 
  DL.concat (DL.concat
            (DL.map (\(_,_,c) -> [[DL.intercalate ":" allvariantinfo] DL.++
                                 [DL.intercalate ":" allregioninfo]  DL.++
                                 [[c]]])
  xs))
    where
      allvariantinfo = DL.map (DText.unpack)
                              (DL.map (vsample)
                              (DL.map (\(x,_,_) -> x) xs))     DL.++
                       DL.map (DText.unpack)
                              (DL.map (vsymbol)
                              (DL.map (\(x,_,_) -> x) xs))     DL.++
                       DL.map (DText.unpack)
                              (DL.map (vchromosome)
                              (DL.map (\(x,_,_) -> x) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (vstartpos)
                              (DL.map (\(x,_,_) -> x) xs))   DL.++
                       DL.map (DText.unpack)
                              (DL.map (vendpos)
                              (DL.map (\(x,_,_) -> x) xs))     DL.++
                       DL.map (DText.unpack)
                              (DL.map (vref)
                              (DL.map (\(x,_,_) -> x) xs))        DL.++
                       DL.map (DText.unpack)
                              (DL.map (valt)
                              (DL.map (\(x,_,_) -> x) xs))        DL.++
                       DL.map (DText.unpack)
                              (DL.map (venst)
                              (DL.map (\(x,_,_) -> x) xs))
      allregioninfo  = DL.map (DText.unpack)
                              (DL.map (rchromosome)
                              (DL.map (\(_,y,_) -> y) xs)) DL.++
                       DL.map (DText.unpack)
                              (DL.map (rtss)
                              (DL.map (\(_,y,_) -> y) xs))        DL.++
                       DL.map (DText.unpack)
                              (DL.map (rstrand)
                              (DL.map (\(_,y,_) -> y) xs))     DL.++
                       DL.map (DText.unpack)
                              (DL.map (rgenename)
                              (DL.map (\(_,y,_) -> y) xs))

--prepareWithinTSS -> This function will
--prepare withintss (see main) for
--printFile function.
prepareWithinTSS :: [[(Variants,BioMartRegion,Char)]] -> [[String]]
prepareWithinTSS []     = []
prepareWithinTSS (x:xs) =
  [convertToListWithinTSS x]
  DL.++ (prepareWithinTSS xs)

{-------------------}


{-Ambiguity Code Functions.-}

--reverseComplementNucleotide -> This function will
--return the reverse complement of a DB.ByteString.
reverseComplementNucleotide :: DB.ByteString -> DB.ByteString
reverseComplementNucleotide currentsequence =
  DB.pack
  (CBUTF8.encode (DL.map (snd)
                  (DL.concatMap
                  (\y -> DL.filter (\(r,_) -> r == y) revcomplementmapping)
                  (CBUTF8.decode
                  (DB.unpack (DB.reverse currentsequence))))))
    where
      revcomplementmapping = [('A','T')
                             ,('T','A')
                             ,('G','C')
                             ,('C','G')]

--smallGrabFastaSequence -> This function will
--grab the correct fasta sequence
--using chromosome information
--in the region file.
smallGrabFastaSequence :: BioMartRegion -> [EDSS.Sequence] -> Vector Char
smallGrabFastaSequence _             []                                = DVU.empty
smallGrabFastaSequence currentregion (currentsequence:restofsequences) =
  if | decodedcurrentsequencechr == ("chr" DL.++ currentregionchr)
     -> DVU.map (EDAC.toChar) currentsequencechars
     | otherwise
     -> smallGrabFastaSequence currentregion restofsequences
    where
      currentsequencechars      = EDSS.characters currentsequence
      decodedcurrentsequencechr = CBUTF8.decode
                                  (DBL.unpack
                                  (EDSS.name currentsequence))
      currentregionchr          = DText.unpack $
                                  rchromosome currentregion

--grabFastaSequence -> This function will
--grab the correct fasta sequence
--using chromosome information
--in the region file.
grabFastaSequence :: FRIConfig -> BioMartRegion -> IO (Either () (Vector Char))
grabFastaSequence config currentregion = do
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" DL.++ (show currenttandd) DL.++ "] "
                         DL.++ "Extracting fasta sequence associated with chromosome: "
                         DL.++ currentregionchr
                         DL.++ ", tss: "
                         DL.++ currentregiontss
                         DL.++ ", strand: "
                         DL.++ currentregionstrand
                         DL.++ ", gene: "
                         DL.++ currentregiongenename
                         DL.++ " ...")
  --Read the file into a strict ByteString.
  inputfastafile <- DB.readFile (DText.unpack $
                                YamlParser.fasta config)
  --Parse fastafile into a [Sequence].
  let pfastafile = parseOnly (EISF.fasta EDAA.DNAI <* endOfInput)
                             inputfastafile
  --Walk through cases of pfastafile.
  case pfastafile of
    Left  _      -> return $ Left () 
    Right cfasta -> return $ Right (smallGrabFastaSequence currentregion
                                                           cfasta) 
    where
      currentregionchr      = DText.unpack (rchromosome currentregion)
      currentregiontss      = DText.unpack (rtss currentregion)
      currentregionstrand   = DText.unpack (rstrand currentregion)
      currentregiongenename = DText.unpack (rgenename currentregion)

--grabRegionSequence -> This function will
--grab the region of the correct chr
--returned from grabFastaSequence.
grabRegionSequence :: Vector Char -> BioMartRegion -> Int -> DB.ByteString
grabRegionSequence currentsequence currentregion currenttsswinsize =
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
      currentregiontss    = DText.unpack (rtss currentregion)
      currentregionstrand = DText.unpack (rstrand currentregion)

--subStrLocationsSmallForward -> This function will
--find the locations for all given substrings
--found using allStrGeneration.
subStrLocationsSmallForward :: FRIConfig -> [String] -> BioMartRegion -> Vector Char -> IO [[Int]]
subStrLocationsSmallForward _      []                                        _             _              = return []
subStrLocationsSmallForward config (currentmappedambstr:restofmappedambstrs) currentregion finalfastafile = do
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" DL.++ (show currenttandd) DL.++ "] "
                         DL.++ "Processing mapped ambiguity code "
                         DL.++ currentmappedambstr
                         DL.++ " ...")
  if | DMaybe.isJust tsswinsize
     -> do _ <- return $ [DBSDFA.indices (DBC.pack currentmappedambstr)
                                         (grabRegionSequence
                                         finalfastafile
                                         currentregion
                                         (read (DText.unpack $ DMaybe.fromJust tsswinsize) :: Int))]
           subStrLocationsSmallForward config
                                       restofmappedambstrs
                                       currentregion
                                       finalfastafile
     | otherwise
     -> do _ <- return $ [DBSDFA.indices (DBC.pack currentmappedambstr)
                                         (grabRegionSequence
                                         finalfastafile
                                         currentregion
                                         2000)]
           subStrLocationsSmallForward config
                                       restofmappedambstrs
                                       currentregion
                                       finalfastafile
    where
      tsswinsize       = tsswindowsize config 

--subStrLocationsSmallReverse -> This function will
--find the locations for all given substrings
--found using allStrGeneration.
subStrLocationsSmallReverse :: FRIConfig -> [String] -> BioMartRegion -> Vector Char -> IO [[Int]]
subStrLocationsSmallReverse _      []                                        _             _              = return []
subStrLocationsSmallReverse config (currentmappedambstr:restofmappedambstrs) currentregion finalfastafile = do
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" DL.++ (show currenttandd) DL.++ "] "
                         DL.++ "Processing mapped ambiguity code " 
                         DL.++ currentmappedambstr 
                         DL.++ " ...")
  if | DMaybe.isJust tsswinsize
     -> do _ <- return $ [DL.map (\a -> (DBC.length
                                        ((grabRegionSequence finalfastafile
                                                             currentregion
                                                             (read (DText.unpack $ DMaybe.fromJust tsswinsize) :: Int)))) - a - 1)
                                        (DBSDFA.indices (DBC.pack currentmappedambstr)
                                        (reverseComplementNucleotide
                                        (grabRegionSequence finalfastafile
                                                            currentregion
                                                            (read (DText.unpack $ DMaybe.fromJust tsswinsize) :: Int))))]
           subStrLocationsSmallReverse config
                                       restofmappedambstrs
                                       currentregion
                                       finalfastafile
     | otherwise
     -> do _ <- return $ [DL.map (\a -> (DBC.length
                                        ((grabRegionSequence finalfastafile
                                                             currentregion
                                                             2000))) - a - 1)
                                        (DBSDFA.indices (DBC.pack currentmappedambstr)
                                        (reverseComplementNucleotide
                                        (grabRegionSequence finalfastafile
                                                            currentregion
                                                            2000)))]
           subStrLocationsSmallReverse config
                                       restofmappedambstrs
                                       currentregion
                                       finalfastafile
    where
      tsswinsize = tsswindowsize config 

--subStrLocations -> This function will
--find the locations for all given substrings
--found using allStrGeneration.
subStrLocations :: FRIConfig -> [String] -> BioMartRegion -> Vector Char -> IO [[Int]]
subStrLocations _      []                     _             _        = return []
subStrLocations config allmappedambiguitystrs currentregion fastaseq = do
  if | DMaybe.isJust tsswinsize
     -> if | currentregionstrand == "-1"
           -> do reversesubstrlocs <- subStrLocationsSmallReverse config
                                                                  allmappedambiguitystrs
                                                                  currentregion
                                                                  fastaseq
                 return $ ((DL.map (DL.map (\i ->
                            ((((read currentregiontss) - (read (DText.unpack $ DMaybe.fromJust tsswinsize) :: Int)) + i) + 2)))
                          reversesubstrlocs)
                          `CPS.using` (CPS.parList CPS.rdeepseq))
           | otherwise
           -> do forwardsubstrlocs <- subStrLocationsSmallForward config
                                                                  allmappedambiguitystrs
                                                                  currentregion
                                                                  fastaseq
                 return $ ((DL.map (DL.map (\i ->
                            ((read currentregiontss) + i)))
                          forwardsubstrlocs)
                          `CPS.using` (CPS.parList CPS.rdeepseq))
     | otherwise
     -> if | currentregionstrand == "-1"
           -> do reversesubstrlocs <- subStrLocationsSmallReverse config
                                                                  allmappedambiguitystrs
                                                                  currentregion
                                                                  fastaseq
                 return $ ((DL.map (DL.map (\i ->
                            ((((read currentregiontss) - 2000) + i) + 2)))
                          reversesubstrlocs)
                          `CPS.using` (CPS.parList CPS.rdeepseq))
           | otherwise
           -> do forwardsubstrlocs <- subStrLocationsSmallForward config
                                                                  allmappedambiguitystrs
                                                                  currentregion
                                                                  fastaseq
                 return $ ((DL.map (DL.map (\i ->
                            ((read currentregiontss) + i)))
                          forwardsubstrlocs)
                          `CPS.using` (CPS.parList CPS.rdeepseq))
    where
      currentregionstrand = DText.unpack $ rstrand currentregion
      currentregiontss    = DText.unpack $ rtss currentregion
      tsswinsize          = tsswindowsize config

--ambiguityCodesWithinRegionCheckSmall -> This function will
--check to see if the ambiguity codes are within genes
--2 kb (or custom TSS window size if provided)
--of the genes TSS.
ambiguityCodesWithinRegionCheckSmall :: FRIConfig -> (String,String) -> [(String,String)] -> [BioMartRegion] -> IO [(String,[String],[String],[[Int]])]
ambiguityCodesWithinRegionCheckSmall _      ([],[])                                        []                     _                             = return []
ambiguityCodesWithinRegionCheckSmall _      _                                              []                     _                             = return []
ambiguityCodesWithinRegionCheckSmall _      ([],[])                                        _                      _                             = return []
ambiguityCodesWithinRegionCheckSmall _      ((_:_),[])                                     (_:_)                  []                            = return []
ambiguityCodesWithinRegionCheckSmall _      ((_:_),(_:_))                                  (_:_)                  []                            = return []
ambiguityCodesWithinRegionCheckSmall _      ([],(_:_))                                     (_:_)                  []                            = return []
ambiguityCodesWithinRegionCheckSmall config currentambtuple@(currentambcode,currentstrand) allmappedambiguitystrs (currentregion:restofregions) = do
  currenttandd <- DTime.getZonedTime
  _ <- SIO.putStrLn ("[" DL.++ (show currenttandd) DL.++ "] "
                         DL.++ "Processing region data associated with gene "
                         DL.++ currentregiongenename
                         DL.++ " ...")  
  if | (ignorestrandedness config) == False
     -> if | currentregionstrand == currentstrand
           -> do --Grab locations of mapped am codes,
                 --and recurse.
                 fastaseq <- grabFastaSequence config
                                               currentregion
                 if | DE.isRight fastaseq
                    -> do let Right finalfastafile = fastaseq
                          substrlocs <- subStrLocations config
                                                        (DL.map (fst) allmappedambiguitystrs)
                                                        currentregion
                                                        finalfastafile                  
                          _ <- return $ [( currentambcode
                                         , DL.map (fst) allmappedambiguitystrs
                                         , allcurrentregiondata
                                         , substrlocs
                                        )]
                          ambiguityCodesWithinRegionCheckSmall config
                                                               currentambtuple
                                                               allmappedambiguitystrs
                                                               restofregions
                    | otherwise
                    -> do _ <- return $ [( currentambcode
                                         , DL.map (fst) allmappedambiguitystrs
                                         , allcurrentregiondata
                                         , [[]]
                                        )]
                          ambiguityCodesWithinRegionCheckSmall config
                                                               currentambtuple
                                                               allmappedambiguitystrs
                                                               restofregions
           | otherwise
           -> --Current ambiguity codes and mapped strings
              --are not correct for region strand
              --(i.e. "-1" != "1" or "1" != "-1").
              ambiguityCodesWithinRegionCheckSmall config
                                                   currentambtuple
                                                   allmappedambiguitystrs
                                                   restofregions
     | otherwise
     -> do --Ignore strandedness, find both the ambiguity mapped strings
           --and the reverse complement ambiguity mapped strings.
           fastaseq <- grabFastaSequence config
                                         currentregion
           if | DE.isRight fastaseq
              -> do let Right finalfastafile = fastaseq
                    substrlocs <- subStrLocations config
                                                  (DL.map (fst) allmappedambiguitystrs)
                                                  currentregion
                                                  finalfastafile
                    _ <- return $ [( currentambcode
                                  , DL.map (fst) allmappedambiguitystrs
                                  , allcurrentregiondata
                                  , substrlocs
                                  )]
                    ambiguityCodesWithinRegionCheckSmall config
                                                         currentambtuple
                                                         allmappedambiguitystrs
                                                         restofregions
              | otherwise
              -> do _ <- return $ [( currentambcode
                                  , DL.map (fst) allmappedambiguitystrs
                                  , allcurrentregiondata
                                  , [[]]
                                  )]
                    ambiguityCodesWithinRegionCheckSmall config
                                                         currentambtuple
                                                         allmappedambiguitystrs
                                                         restofregions
    where
      allcurrentregiondata  = [ currentregionchr
                              , currentregiontss
                              , currentregionstrand
                              , currentregiongenename
                              ] 
      currentregionchr      = DText.unpack (rchromosome currentregion)
      currentregiontss      = DText.unpack (rtss currentregion)
      currentregionstrand   = DText.unpack (rstrand currentregion)
      currentregiongenename = DText.unpack (rgenename currentregion)

--ambiguityCodesWithinRegionCheck -> This function will
--check to see if the ambiguity codes are within genes
--2 kb (or custom TSS window size if provided)
--of the genes TSS.
ambiguityCodesWithinRegionCheck :: FRIConfig -> [(String,String)] -> [[(String,String)]] -> [BioMartRegion] -> IO [[(String,[String],[String],[[Int]])]]
ambiguityCodesWithinRegionCheck _      []                              []                                                              _          = return []
ambiguityCodesWithinRegionCheck _      _                               []                                                              _          = return []
ambiguityCodesWithinRegionCheck _      []                              _                                                               _          = return []
ambiguityCodesWithinRegionCheck config (currentambcode:restofambcodes) (currentmappedambiguitystrgroup:restofmappedambiguitystrgroups) allregions = do
  ambcodeswithinregioncheck <- ambiguityCodesWithinRegionCheckSmall config
                                                                    currentambcode
                                                                    currentmappedambiguitystrgroup
                                                                    allregions
  _ <- return $ [ambcodeswithinregioncheck]
  ambiguityCodesWithinRegionCheck config
                                  restofambcodes
                                  restofmappedambiguitystrgroups
                                  allregions
 
--allStrGenerationSmall -> This function will
--return all substring locations with current
--ambiguity code mapped to true nucleotides.
allStrGenerationSmall :: String -> IO [String]
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
      nucleotidemapping = [('A','A')
                          ,('T','T')
                          ,('G','G')
                          ,('C','C')
                          ,('Y','C')
                          ,('Y','T')
                          ,('R','A')
                          ,('R','G')
                          ,('W','A')
                          ,('W','T')
                          ,('S','G')
                          ,('S','C')
                          ,('K','T')
                          ,('K','G')
                          ,('M','C')
                          ,('M','A')
                          ,('D','A')
                          ,('D','G')
                          ,('D','T')
                          ,('V','A')
                          ,('V','C')
                          ,('V','G')
                          ,('H','A')
                          ,('H','C')
                          ,('H','T')
                          ,('B','C')
                          ,('B','G')
                          ,('B','T')
                          ,('N','A')
                          ,('N','T')
                          ,('N','G')
                          ,('N','C')
                          ,('X','A')
                          ,('X','T')
                          ,('X','G')
                          ,('X','C')] 

--allStrGeneration -> This function will
--return all substring locations with all
--ambiguity code mapped to true nucleotides.
allStrGeneration :: [String] -> IO [[String]]
allStrGeneration [] = return []
allStrGeneration xs = DT.mapM (allStrGenerationSmall) xs

--ambiguityCodesReverseComplementSmall -> This function will
--calculate the reverse complement for user specified 
--ambiguity codes.
ambiguityCodesReverseComplementSmall :: [String] -> [String]
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
      revcomplementmapping = [('A','T')
                             ,('T','A')
                             ,('G','C')
                             ,('C','G')
                             ,('Y','R')
                             ,('R','Y')
                             ,('W','W')
                             ,('S','S')
                             ,('K','M')
                             ,('M','K')
                             ,('D','H')
                             ,('H','D')
                             ,('V','B')
                             ,('B','V')
                             ,('X','X')
                             ,('N','N')
                             ,('-','-')]  

--ambiguityCodesReverseComplement -> This function will
--grab an ambiguity codes reverse complement.
ambiguityCodesReverseComplement :: FRIConfig -> [String]
ambiguityCodesReverseComplement config =
  ambiguityCodesReverseComplementSmall ambcodes  
    where
      ambcodes = DL.map (DText.unpack)
                 (ambiguitycodes config)

--tupleConverterAmbCodesWithinTSS -> This function will
--convert a 4-tuple to the corrected
--4-tuple.
tupleConverterAmbCodesWithinTSS :: (String,[String],[String],[[Int]]) -> ([[String]],[[String]],[[String]],[[Int]])
tupleConverterAmbCodesWithinTSS (a,b,c,d) =
  ((DL.map (\x -> [x])
   (DL.replicate (DL.length b) a))
  ,(DL.map (\x -> [x]) b)
  ,(DL.replicate (DL.length b) c)
  ,d)

--tupleToListAmbCodesWithinTSS -> This function will
--convert a 4-tuple to a list.
tupleToListAmbCodesWithinTSS :: ([[String]],[[String]],[[String]],[[Int]]) -> [[String]]
tupleToListAmbCodesWithinTSS ([],    _,     _,     _)      = []
tupleToListAmbCodesWithinTSS ((_:_), [],    _,     _)      = []
tupleToListAmbCodesWithinTSS ((_:_), (_:_), [],    _)      = []
tupleToListAmbCodesWithinTSS ((_:_), (_:_), (_:_), [])     = []
tupleToListAmbCodesWithinTSS ((a:as),(b:bs),(c:cs),(d:ds)) =
  [a DL.++ b DL.++ c DL.++ (DL.map (show) d)]
  DL.++ (tupleToListAmbCodesWithinTSS (as,bs,cs,ds))

--convertToListAmbCodesWithinTSS -> This function will
--convert 4-tuple to list.
convertToListAmbCodesWithinTSS :: [(String,[String],[String],[[Int]])] -> [[[String]]]
convertToListAmbCodesWithinTSS [] = []
convertToListAmbCodesWithinTSS xs =
  DL.map (\y -> tupleToListAmbCodesWithinTSS y)
         (DL.map (\x -> tupleConverterAmbCodesWithinTSS x)
  xs)

--prepareAmbiguityCodesWithinTSS -> This function will
--prepare ambiguitycodeswithintss (see main) for
--printFile function.
prepareAmbiguityCodesWithinTSS :: [[(String,[String],[String],[[Int]])]] -> [[[String]]]
prepareAmbiguityCodesWithinTSS []     = []
prepareAmbiguityCodesWithinTSS (x:xs) =
  (convertToListAmbCodesWithinTSS x)
  DL.++ (prepareAmbiguityCodesWithinTSS xs)

{---------------------------}


{-Variant functions.-}

--variantsAmbiguityCodesCheckerSmaller -> This function will
--check whether the variant in question lies within an
--ambiguity code sequence.
variantsAmbiguityCodesCheckerSmaller :: (Variants,BioMartRegion,Char) -> [String] -> Int -> [String]
variantsAmbiguityCodesCheckerSmaller _                                   []     _ = []
variantsAmbiguityCodesCheckerSmaller xs@(currentvariant,currentregion,_) (y:ys) z =
  --TSS reads in reverse direction (-1).
  if | currentregionstrand == "-1"
     -> if | (read y :: Int) >= (read currentvariantstartpos :: Int) &&
             (read currentvariantstartpos :: Int) >= ((((read y) - z) + 1) :: Int)
           ->[y] DL.++ (variantsAmbiguityCodesCheckerSmaller xs ys z)
           | otherwise
           -> variantsAmbiguityCodesCheckerSmaller xs ys z
     --TSS reads in forward direction (1).
     | (read y :: Int) <= (read currentvariantstartpos :: Int) &&
       (read currentvariantstartpos :: Int) <= ((((read y) + z) - 1) :: Int)
     -> [y] DL.++ (variantsAmbiguityCodesCheckerSmaller xs ys z)
     | otherwise
     -> variantsAmbiguityCodesCheckerSmaller xs ys z
    where
      currentregionstrand    = DText.unpack $ rstrand currentregion
      currentvariantstartpos = DText.unpack $ vstartpos currentvariant      

--variantsAmbiguityCodesCheckerSmall -> This function will
--call variantsAmbiguityCodesCheckerSmaller.
variantsAmbiguityCodesCheckerSmall :: (Variants,BioMartRegion,Char) -> [[String]] -> [[String]]
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
      allregiondata  = [ DText.unpack $ rchromosome currentregion
                       , DText.unpack $ rtss currentregion
                       , DText.unpack $ rstrand currentregion
                       , DText.unpack $ rgenename currentregion
                       ] 
      allvariantdata = [ DText.unpack $ vsample currentvariant
                       , DText.unpack $ vsymbol currentvariant
                       , DText.unpack $ vchromosome currentvariant
                       , DText.unpack $ vstartpos currentvariant
                       , DText.unpack $ vendpos currentvariant
                       , DText.unpack $ vref currentvariant
                       , DText.unpack $ valt currentvariant
                       , DText.unpack $ venst currentvariant
                       ]

--variantsAmbiguityCodesChecker -> This function will
--call variantsAmbiguityCodesCheckerSmall.
variantsAmbiguityCodesChecker :: (Variants,BioMartRegion,Char) -> [[[String]]] -> [[String]]
variantsAmbiguityCodesChecker _      []     = []
variantsAmbiguityCodesChecker xs     (y:ys) =
  (variantsAmbiguityCodesCheckerSmall xs y)
  DL.++ (variantsAmbiguityCodesChecker xs ys)

--variantsWithinAmbiguityCodesAndTSSSmall -> This function will
--grab all filtered ambiguity codes for matching genes.
variantsWithinAmbiguityCodesAndTSSSmall :: [(Variants,BioMartRegion,Char)] -> [[[String]]] -> [[String]]
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
      currentvariantsymbol = DText.unpack $ vsymbol currentvariant

--variantsWithinAmbiguityCodesAndTSS ->  This function will
--identify all variants that are within ambiguity codes
--and corresponding genes TSS.
variantsWithinAmbiguityCodesAndTSS :: [[(Variants,BioMartRegion,Char)]] -> [[[String]]] -> [[String]]
variantsWithinAmbiguityCodesAndTSS []     [] = []
variantsWithinAmbiguityCodesAndTSS _      [] = []
variantsWithinAmbiguityCodesAndTSS []     _  = []
variantsWithinAmbiguityCodesAndTSS (x:xs) ys =
  (variantsWithinAmbiguityCodesAndTSSSmall variantsfiltered ys) 
  DL.++ (variantsWithinAmbiguityCodesAndTSS xs ys)
    where 
      variantsfiltered = DL.filter (\(_,_,c) -> c == 'Y') x 

{--------------------}


{-To CSV function.-}

toCSV :: [[String]] -> String
toCSV allrows = genCsvFile allrows

{------------------}


{-CSV printing functions.-}

writeCSV :: FRIConfig -> String -> String -> IO ()
writeCSV _      []      []             = return ()
writeCSV config csvdata outputfilename =
  if | DL.last outputdir == '/'
     -> SIO.writeFile (outputdir DL.++
                       outputfilename)
                       csvdata
     | otherwise
     -> SIO.writeFile (outputdir DL.++
                      "/"        DL.++ 
                      outputfilename)
                      csvdata
    where
      outputdir = DText.unpack $ outputdirectory config      

{-------------------------}
