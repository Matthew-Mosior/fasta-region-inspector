{-# LANGUAGE DataKinds                 #-}
{-# LANGUAGE FlexibleContexts          #-}
{-# LANGUAGE LinearTypes               #-}
{-# LANGUAGE MultiWayIf                #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE OverloadedStrings         #-}
{-# LANGUAGE RankNTypes                #-}
{-# LANGUAGE TypeOperators             #-}
{-# LANGUAGE TypeFamilies              #-}

module DiagramGeneration.Generate where

import Logging
import Types
import YamlParser

import Linear.UtilityLinear

import           Codec.Binary.UTF8.String        as CBUTF8
import qualified Control.Lens                    as CL
import           Data.ByteString                 as DB
import           Data.List                       as DL
import           Data.List.Split                 as DLS
import           Data.Maybe
import           Data.Text                       as DText
import           Data.Text.Encoding              as DTE
import           Data.Traversable                as DT
import           Diagrams.Combinators
import           Diagrams.Prelude
import           Diagrams.TwoD
import           Diagrams.TwoD.Path
import           Diagrams.TwoD.Text
import           Diagrams.TwoD.Types
import           Data.Typeable
import           Diagrams.Size
import           Diagrams.Backend.SVG
import           Effectful
import           Effectful.Ki
import           GHC.Stack
import           Graphics.SVGFonts
import qualified System.IO.Resource.Linear       as Linear

generateGraphics :: forall {es :: [Effect]} {b}.
                    ( FRILogging :> es
                    , StructuredConcurrency :> es
                    , IOE :> es
                    , HasCallStack
                    )
                 => FRIConfig
                 -> [[String]]
                 -> Eff es ()
generateGraphics config
                 finalvariants = do
  let ((callingfunction,_):_) = getCallStack callStack
  _                           <- showPrettyLog LogInfo
                                               (maxnumberconcthreads config)
                                               callingfunction
                                               "Filtering out variants that are found \
                                               \within a mapped ambiguity string, \
                                               \and that also lie within 2 kb of the TSS from the final data."
  let finalvariants'          = DL.filter (\xs -> xs DL.!! 3 /= "N/A" &&
                                                  xs DL.!! 4 /= "N/A" &&
                                                  xs DL.!! 5 /= "N/A"
                                          )
                                finalvariants
  scoped $ \scope -> do
    DT.forM finalvariants' $ \currentfinalvariant -> do
      fork scope (do let currentregion = BioMartRegion { biomartregion_sequencedescription = DText.pack (DLS.splitOn ":" (currentfinalvariant DL.!! 1) DL.!! 0)
                                                       , biomartregion_tss                 = DText.pack (DLS.splitOn ":" (currentfinalvariant DL.!! 1) DL.!! 1)
                                                       , biomartregion_strand              = DText.pack (DLS.splitOn ":" (currentfinalvariant DL.!! 1) DL.!! 2)
                                                       , biomartregion_genename            = DText.pack (DLS.splitOn ":" (currentfinalvariant DL.!! 1) DL.!! 3)
                                                       }
                     let currentregiongenename = DText.unpack (biomartregion_genename currentregion)
                     _                 <- showPrettyLog LogInfo
                                                        (maxnumberconcthreads config)
                                                        callingfunction
                                                        ( "Reading fasta index file for "
                                                          DL.++
                                                          currentregiongenename
                                                          DL.++
                                                          "."
                                                        )
                     cfai              <- liftIO $ Linear.run $ getFAILineLinear config
                                                                                 currentregion
                     case cfai of
                       Nothing -> do _ <- showPrettyLog LogTrace
                                                        (maxnumberconcthreads config)
                                                        callingfunction
                                                        ( "Could not retrieve required information from fasta index file for "
                                                          DL.++
                                                          currentregiongenename
                                                          DL.++
                                                          "."
                                                        )
                                     return ()
                       Just cfaif -> do _ <- showPrettyLog LogInfo
                                                           (maxnumberconcthreads config)
                                                           callingfunction
                                                           ( "Reading in fasta file for "
                                                             DL.++
                                                             currentregiongenename
                                                             DL.++
                                                             "."
                                                           )
                                        let ambiguitycodelocation = currentfinalvariant DL.!! 5
                                        fastaseq   <- liftIO $ Linear.run $ getFASTASequenceGraphicLinear config
                                                                                                          (read ambiguitycodelocation :: Int)
                                                                                                          cfaif
                                        case fastaseq of
                                          Nothing        -> do _ <- showPrettyLog LogTrace
                                                                                  (maxnumberconcthreads config)
                                                                                  callingfunction
                                                                                  ( "Could not retrieve required information from fasta file for "
                                                                                    DL.++
                                                                                    currentregiongenename
                                                                                    DL.++
                                                                                    "."
                                                                                    )
                                                               return ()
                                          Just fastaseqf -> do let graphicfilename   = (DLS.splitOn ":" (currentfinalvariant DL.!! 0) DL.!! 0) DL.++
                                                                                       "_"                                                     DL.++
                                                                                       (DLS.splitOn ":" (currentfinalvariant DL.!! 0) DL.!! 2) DL.++
                                                                                       ":"                                                     DL.++
                                                                                       (DLS.splitOn ":" (currentfinalvariant DL.!! 0) DL.!! 3) DL.++
                                                                                       "-"                                                     DL.++
                                                                                       (DLS.splitOn ":" (currentfinalvariant DL.!! 0) DL.!! 4) DL.++
                                                                                       "_"                                                     DL.++
                                                                                       (DLS.splitOn ":" (currentfinalvariant DL.!! 0) DL.!! 5) DL.++
                                                                                       ":"                                                     DL.++
                                                                                       (DLS.splitOn ":" (currentfinalvariant DL.!! 0) DL.!! 6)
                                                               _                     <- showPrettyLog LogInfo
                                                                                                      (maxnumberconcthreads config)
                                                                                                      callingfunction
                                                                                                      ( "Generating graphic for "
                                                                                                        DL.++
                                                                                                        graphicfilename
                                                                                                        DL.++
                                                                                                        "."
                                                                                                      )
                                                               let sequencetoprint   = annotatedFASTASequence fastaseqf
                                                                                                              currentfinalvariant
                                                                                                              currentregion
                                                               let sequencetoprintf  = DL.map (\(a,b) -> case (biomartregion_strand $ currentregion) of
                                                                                                           "1" -> do let stp = DText.unpack     $
                                                                                                                                 decodeUtf8     $
                                                                                                                                   DB.singleton b
                                                                                                                     if | a == Just VariantIndex
                                                                                                                        -> text stp                  #
                                                                                                                           fc ambiguitycodecolor     #
                                                                                                                           fontSizeL variantfontsize #
                                                                                                                           ultraBold                 #
                                                                                                                           oblique
                                                                                                                        | a == Just AmbiguityStringIndice
                                                                                                                        -> text stp                        #
                                                                                                                           fc ambiguitycodecolor           #
                                                                                                                           fontSizeL ambiguitycodefontsize #
                                                                                                                           ultraBold                       #
                                                                                                                           oblique 
                                                                                                                        | isNothing a &&
                                                                                                                          stp == "A"
                                                                                                                        -> text stp        #
                                                                                                                           fc adeninecolor #
                                                                                                                           fontSizeL nonambiguitycodefontsize
                                                                                                                        | isNothing a &&
                                                                                                                          stp == "T"
                                                                                                                        -> text stp        #
                                                                                                                           fc thyminecolor #
                                                                                                                           fontSizeL nonambiguitycodefontsize
                                                                                                                        | isNothing a &&
                                                                                                                          stp == "G"
                                                                                                                        -> text stp        #
                                                                                                                           fc guaninecolor #
                                                                                                                           fontSizeL nonambiguitycodefontsize
                                                                                                                        | otherwise
                                                                                                                        -> text stp         #
                                                                                                                           fc cytosinecolor #
                                                                                                                           fontSizeL nonambiguitycodefontsize
                                                                                                           _   -> do let stp = DText.unpack                    $
                                                                                                                                 decodeUtf8                    $
                                                                                                                                   reverseComplementNucleotide $
                                                                                                                                     DB.singleton b
                                                                                                                     if | a == Just VariantIndex
                                                                                                                        -> text stp                  #
                                                                                                                           fc ambiguitycodecolor     #
                                                                                                                           fontSizeL variantfontsize #
                                                                                                                           ultraBold                 #
                                                                                                                           oblique
                                                                                                                        | a == Just AmbiguityStringIndice
                                                                                                                        -> text stp                        #
                                                                                                                           fc ambiguitycodecolor           #
                                                                                                                           fontSizeL ambiguitycodefontsize #
                                                                                                                           ultraBold                       #
                                                                                                                           oblique 
                                                                                                                        | isNothing a &&
                                                                                                                          stp == "A"
                                                                                                                        -> text stp        #
                                                                                                                           fc adeninecolor #
                                                                                                                           fontSizeL nonambiguitycodefontsize
                                                                                                                        | isNothing a &&
                                                                                                                          stp == "T"
                                                                                                                        -> text stp        #
                                                                                                                           fc thyminecolor #
                                                                                                                           fontSizeL nonambiguitycodefontsize
                                                                                                                        | isNothing a &&
                                                                                                                          stp == "G"
                                                                                                                        -> text stp        #
                                                                                                                           fc guaninecolor #
                                                                                                                           fontSizeL nonambiguitycodefontsize
                                                                                                                        | otherwise
                                                                                                                        -> text stp         #
                                                                                                                           fc cytosinecolor #
                                                                                                                           fontSizeL nonambiguitycodefontsize
                                                                                              )
                                                                                       ((\(FASTASequenceGraphics fsg) -> fsg) sequencetoprint)
                                                               let currentgraphic    = sequencetoprintf            #
                                                                                       hsep 38                     #
                                                                                       fc black                    #
                                                                                       centerXY
                                                               currentgraphicf       <- return $
                                                                                          currentgraphic
                                                                                          `atop` 
                                                                                          ( rect (fromIntegral 800) (fromIntegral 400) #
                                                                                            bg white
                                                                                          )
                                                               liftIO $ renderSVG ( (DText.unpack $ outputdirectory config) DL.++
                                                                                    graphicfilename                         DL.++
                                                                                    ".svg"
                                                                                  )
                                                                                  (dims $ V2 800 400)
                                                                                  currentgraphicf
                 )
    atomically $ awaitAll scope
  where
    variantfontsize          = 44.5 :: Double
    ambiguitycodefontsize    = 38.5 :: Double
    nonambiguitycodefontsize = 32.5 :: Double
    titlefontsize            = 30.5 :: Double
    adeninecolor             = limegreen
    thyminecolor             = firebrick
    guaninecolor             = orangered 
    cytosinecolor            = midnightblue
    ambiguitycodecolor       = black 
    annotatedFASTASequence :: FASTASequence
                           -> [String]
                           -> BioMartRegion
                           -> FASTASequenceGraphics
    annotatedFASTASequence fs
                           currentfinalvariant
                           currentregion= do
      let ambiguitycodelength   = DL.length $
                                    currentfinalvariant DL.!! 4
      let ambiguitycodeelements = case (biomartregion_strand $ currentregion) of
                                    "1" -> [9..((9+ambiguitycodelength)-1)]
                                    _   -> [(9-ambiguitycodelength+1)..9]
      let variantstring         = ( DLS.splitOn ":" $
                                    currentfinalvariant DL.!! 0
                                  ) DL.!! 3
      let variantelement        = do let variantelement' = (read variantstring) - (read $ currentfinalvariant DL.!! 5)
                                     if | variantelement' < 0
                                        -> [9-(abs variantelement')]
                                        | otherwise
                                        -> [9+variantelement']
      let fs'                   = DB.unpack $ (\(FASTASequence fs') -> fs') fs
      let fs''                  = DL.zip [0..(DL.length fs')] fs'
      let fs'''                 = DL.map (\(a,b) -> if | a `DL.elem` ambiguitycodeelements &&
                                                         a `DL.notElem` variantelement
                                                       -> (Just AmbiguityStringIndice,b)
                                                       | a `DL.elem` ambiguitycodeelements &&
                                                         a `DL.elem` variantelement
                                                       -> (Just VariantIndex,b)
                                                       | otherwise
                                                       -> (Nothing,b)
                                         )
                                  fs''
      FASTASequenceGraphics fs'''
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
