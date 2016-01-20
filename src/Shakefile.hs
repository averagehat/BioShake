{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE FlexibleContexts #-} -- needed for Regex (=~)
module Lib
    ( block
    ) where

import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util
import Text.Regex.Posix
import Data.String.Here
import Control.Monad
import Control.Applicative
-- To run:
-- stack setup
-- stack exec shell-exe

block = shakeArgs shakeOptions{shakeFiles="_build"} $ do

    want ["data/align/align.bam"] 

    "data/align/align.bam" %> \out -> do -- align.bam is the target
       fqs <- allFqs
       let trimmed = map (<.> "cutadapt") fqs
       need trimmed
       --TODO: replace with mapper
       cmd Shell "cat" trimmed ">" out 

    "data/*.cutadapt" %> \out -> do 
       sffs <- getDirectoryFiles "" ["data/*.sff"]
       let fqs = [c -<.> "fastq" | c <- sffs]  -- -<.> replaces extension
       need fqs
       allFastqs <- getDirectoryFiles "" ["data/*.fastq"] 
       mapPairedUnpaired pCutAdapt unpCutAdapt  allFastqs 

    "data/*.fastq" %> \out -> do 
         need [out -<.> "sff"] 
         runPython [i|
         from Bio import SeqIO
         SeqIO.convert('${src}', 'sff', '${out}', 'fastq')
         |]

unpCutAdapt fq = unit $ cmd "cutadapt" ["-a", "AACCGGTT", "-o", fq <.> "cutadapt" , fq ] 
pCutAdapt fwd rev = unit $ cmd "cutadapt" ["-a", "ADAPTER_FWD", "-A", "ADAPTER_REV", "-o", outFwd, "-p", outRev,  fwd, rev]
  where (outFwd, outRev) = (fwd <.> "cutdapt", rev <.> "cutadapt")

matching str strings = (filter (=~ str) strings) :: [String] 

groupFastqs fqs = (unpaired, fwd, rev) 
  where
    fwd = matching "_R1_" fqs
    rev = matching "_R2_" fqs 
    unpaired = [x | x <- fqs, not (x `elem` (fwd ++ rev))]

-- FilePath is an alias for string
mapPairedUnpaired :: (FilePath -> FilePath -> Action ()) -> (FilePath -> Action ()) -> [FilePath] -> Action ()
mapPairedUnpaired pF unpF fqs = do
  _ <- (mapM_ unpF unp) 
  zipWithM_ pF fwd rev 
  where (unp, fwd, rev) =  groupFastqs fqs

 
allFqs = liftA2 (++) sffs fqs
  where 
    sffs = liftM (map (-<.> "fastq")) $ getDirectoryFiles "" ["data/*.sff"]
    fqs =  getDirectoryFiles "" ["data/*.fastq"] 

runPython str = do 
  withTempFile $ \file -> do 
    liftIO $ writeFile file fixedStr
    cmd "python" file 
  where 
    isWhitespace = null . dropWhile (\x -> (x==' ') || (x == '\n'))
    full = dropWhile isWhitespace $ lines str 
    indent = length $ takeWhile (== ' ') $ head full
    fixedStr = unlines $ map (drop indent) full

-- note that stack cannot tell if a rule has changed. so altering want/need and re-running won't work.
       -- () <-  is needed if cmd is not the last statement, e.g.
       -- () <- cmd Shell "cat" fqs ">" [out]  -- [out] = ["align.bam"]
         -- sffs <- getDirectoryFiles "" ["//*.sff"]
    --foldl (\x y -> x ++ "\n" ++ y) "" $ map (drop indent) $ lines str
