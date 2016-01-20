{-# LANGUAGE QuasiQuotes #-} 
--{-# LANGUAGE OverlappingInstances #-} --lets override tuple

import Data.Char(ord,chr,toUpper)
import Control.Monad.Trans(liftIO)
import Data.String.Here
import System.Directory(doesFileExist)
import Text.Regex(subRegex,mkRegex)
import Control.Applicative
import Control.Monad.Extra(ifM)

data IUPAC = A | C | G | T | N
  deriving (Eq, Ord, Show, Enum) --[A ..] == [A, C, G, T]

type Id = String
type Quality = [Int] 
-- data FastaRecord = Rec Id DNASeq
type DNASeq = [IUPAC] 
data FastqRecord = Rec Id DNASeq Quality
instance Show FastqRecord  where 
     show (Rec id seq qual) = [i|@${id}
+
${join $ map show seq}
${map encodeQual qual}|] 
       where join = foldr (++) "" 

toIUPAC nt = [A ..] !! (length $ takeWhile (/= (toUpper nt)) "ACGTN")
decodeQual = (subtract 33) . ord
encodeQual = chr . (33 +)
makeSeq (a:b:_:d:[]) = Rec (tail a) (map toIUPAC b) (map decodeQual d)

splitEvery n = takeWhile (not . null) . map (take n) . iterate (drop n)

readFastq :: FilePath -> IO [FastqRecord] 
--e.g. fmap head $ readFastq "foo.fastq"
readFastq p = do
  str <- readFile p
  return $ map makeSeq $ splitEvery 4 $ lines str
  -- <- NB: is like fmap
  --
  --
--  "data/*.filtered" %> \out -> do
--     let src = out -<.> "fastq"
--     need [src]
--     let _id = getIndex src 
--     seqs <- readFastq src
--     filtered <- (dropNs _id) $ filter belowQual seqs
--     -- write filtered to file


getIndex p = do 
  exists <- doesFileExist name
  if exists && (name /= p) then return (Just name) else return Nothing
  where name = subRegex (mkRegex "_R([12])_") p "_I\\1_"

belowQual n (Rec _ _ qual) = (minimum qual) < n 

dropNs :: IO [FastqRecord] -> Maybe FilePath -> IO [FastqRecord]
dropNs recs Nothing = recs
dropNs recs (Just _id) = do
    idSeqs <- (readFastq _id)
    seqs  <- recs
    return $ map snd $ filter (hasN . fst) $ zip idSeqs seqs
 where 
   hasN (Rec _ seq _) = N `elem` seq


-- filter (belowQual n) 
--map snd $ filter (hasN . fst) $ zip (readFastq _id) recs
ifIO p s1 s2 = do 
  b <- p 
  if b then return s1 else return s2
