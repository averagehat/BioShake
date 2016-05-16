{-# LANGUAGE QuasiQuotes #-} 
module Bio where
import Data.Char(ord,chr,toUpper)
import Data.String.Here
import System.Directory(doesFileExist)
import Text.Regex(subRegex,mkRegex)
import Control.Applicative

data IUPAC = A | C | G | T | N
  deriving (Eq, Ord, Show, Enum) --[A ..] == [A, C, G, T]

type Id = String
type Quality = String
type DNASeq = [IUPAC] 
data FastqRecord = Rec Id DNASeq Quality

toStr :: FastqRecord -> String
toStr (Rec id seq qual) = [i|@${id}
+
${foldr1 (++) $ map show seq}
${qual}|] 


toIUPAC :: Char -> IUPAC
toIUPAC nt = [A ..] !! (length $ takeWhile (/= (toUpper nt)) "ACGTN")

--TODO: make parsing safe
makeSeq :: [String] -> FastqRecord
makeSeq (a:b:_:d:[]) = Rec (tail a) (map toIUPAC b) d

readFastq :: FilePath -> IO [FastqRecord] 
readFastq p = do
  str <- readFile p
  return $ map makeSeq $ chunksOf 4 $ lines str
  where chunksOf n = takeWhile (not . null) . map (take n) . iterate (drop n)

-- get the index string for a specific read
-- >>> getIndex "foo_R1_.fastq"
-- foo_I1_.fastq
getIndex :: FilePath -> IO (Maybe FilePath)
getIndex p = do
  exists <- doesFileExist name
  return $ if exists && (name /= p) then (Just name) else Nothing
  where name = subRegex (mkRegex "_R([12])_") p "_I\\1_"

-- drop Ns and reads with bad quality index.
applyFilters :: Int -> FilePath -> IO [FastqRecord]
applyFilters q fp = fmap dropNs (dropIndexQual q fp)

-- filter out reads where the read sequence contains an "N"
dropNs seqs = filter (not . hasN) seqs
  where hasN (Rec _ seq _) = N `elem` seq
  
-- Given a miseq file fp, and Integer N
-- 1. find its matching index file
-- 2. drop all reads from fp where the matching index read has quality under N
dropIndexQual :: Int -> FilePath -> IO [FastqRecord]
dropIndexQual n fp = do
 index  <- getIndex fp
 seqs   <- readFastq fp
 let idSeqs = fmap readFastq index
 maybe (return seqs) (liftA (dropBad seqs)) idSeqs
 where 
  dropBad recs idSeqs' = map snd $ filter (goodQual . fst) $ zip idSeqs' recs
  goodQual (Rec _ _ qual) = (minimum $ map decodeQual qual) >= n
    where decodeQual = (subtract 33) . ord

--encodeQual = chr . (33 +)
--type Quality = [Int] 
