{-# LANGUAGE QuasiQuotes #-} 
module Bio where
import Data.Char(ord,chr,toUpper)
import Data.String.Here
import Data.Word (Word8)
import System.Directory(doesFileExist)
import Text.Regex(subRegex,mkRegex)
import Control.Applicative
import Bio.Sequence.FastQ (readSangerQ, writeSangerQ)
import Bio.Core.Sequence (BioSeqQual, unSD, seqheader, seqdata, seqqual, unQD, unQual)
import qualified Data.ByteString.Lazy.Char8 as C
import qualified Data.ByteString.Lazy as BB 


-- | get the index string for a specific read
-- >>> getIndex "test/data/foo_R1_.fastq"
-- Just "test/data/foo_I1_.fastq" 
getIndex :: FilePath -> IO (Maybe FilePath)
getIndex p = do
  -- TODO: this should only operate on the basename
  exists <- doesFileExist name
  return $ if exists && (name /= p) then (Just name) else Nothing
  where name = subRegex (mkRegex "_R([12])_") p "_I\\1_"

-- drop Ns and reads with bad quality index.
--applyFilters :: BioSeqQual s => Int -> FilePath -> IO [s]
applyFilters q fp = dropNs <$> dropIndexQual q fp
--writeSangerQ "foo" <$> 
-- filter out reads where the read sequence contains an "N"
dropNs :: BioSeqQual s => [s] -> [s]
dropNs seqs = filter (not . hasN) seqs
  where hasN  = ('N' `C.elem`) . unSD . seqdata
  
-- Given a miseq file fp, and Integer N
-- 1. find its matching index file
-- 2. drop all reads from fp where the matching index read has quality under N
--dropIndexQual :: BioSeqQual s => Int -> FilePath -> IO [s]
dropIndexQual n fp = do
 index  <- getIndex fp
 seqs   <- readSangerQ fp
 let idSeqs = fmap readSangerQ index
 maybe (return seqs) (liftA (dropBad seqs)) idSeqs
 where 
  dropBad recs idSeqs' = map snd $ filter (goodQual . fst) $ zip idSeqs' recs
  goodQual seq = (fromIntegral $ BB.minimum $ unQD $ seqqual seq) >= n
