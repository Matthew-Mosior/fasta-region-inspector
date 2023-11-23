module SequenceToVector where

import YamlParser

import Codec.Binary.UTF8.String           as CBUTF8
import Data.ByteString.Lazy               as DBL
import Data.List                          as DL
import Data.Text                          as DText
import Data.Vector.Unboxed                as DVU
import qualified ELynx.Alphabet.Character as EDAC
import qualified ELynx.Sequence.Sequence  as EDSS

extractNameAndCharacters :: [EDSS.Sequence]
                         -> [(Text,Vector Char)]
extractNameAndCharacters allsequences =
  DL.map (\x -> ( DText.pack    $
                  CBUTF8.decode $
                  DBL.unpack    $
                  EDSS.name x
                , DVU.map (EDAC.toChar)
                  (EDSS.characters x)
                )
         )
  allsequences
