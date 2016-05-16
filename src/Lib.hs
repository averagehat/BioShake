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
import Data.Traversable (for)
import Data.Foldable (for_)
import Bio.Sequence.FastQ (writeSangerQ)
import Bio (applyFilters)

-- To run:
-- stack setup
-- stack build
-- stack exec shell-exe


-- -<.> replace the extension
-- <.>  add an extension
-- </>  join paths

block = shakeArgs shakeOptions{shakeFiles="_build"} $ do

    want ["data/align/align.bam"] -- final result of build

    "data/align/align.bam" %> \out -> do -- align.bam is the target, assinged to `out`
      fqs <- allFqs
      let trimmed = map (<.> "cutadapt") fqs
      need trimmed
      --TODO: replace "cat" with mapper
      cmd Shell "cat" trimmed ">" out 

    "data/*.cutadapt" %> \out -> do -- trim with cutadapt.
      fqs <- allFqs
      need fqs
      mapPairedUnpaired pCutAdapt unpCutAdapt fqs

    "data/*.filtered" %> \out -> do
      let src = out -<.> "fastq"
      need [src]
      recs <- liftIO $ applyFilters 30 src
      liftIO $ writeSangerQ out recs
--liftIO $ writeSangerQ out =<< applyFilters 30 src
      
    "data/*.fastq" ?%> \out -> do  -- use python to convert fastq to sff. (?%>) indicates it may already be supplied by user.
      let src = out -<.> "sff"
      need [src]
      -- below is a multiline string
      runPython [i|
      from Bio import SeqIO
      SeqIO.convert('${src}', 'sff', '${out}', 'fastq')
      |]
      
 where      
  -- rule for targets which may already exist.
  -- otherwise, shake will try to build something--even if the user has supplied it.
  (?%>) pat act = pat %> \out -> do
    b <- doesFileExist out
    unless b $ act out
      
-- cutadapt runner functions
unpCutAdapt fq = unit $ cmd "cutadapt" ["-a", "AACCGGTT", "-o", fq <.> "cutadapt" , fq ] 
pCutAdapt fwd rev = unit $ cmd "cutadapt" ["-a", "ADAPTER_FWD", "-A", "ADAPTER_REV", "-o", outFwd, "-p", outRev,  fwd, rev]
  where (outFwd, outRev) = (fwd <.> "cutdapt", rev <.> "cutadapt") 

-- split files based on their names into unpaired, etc.
-- groupFastqs ["unp.fastq", "bar_R1_.fastq", "foo_R1_.fastq", "foo_R2_.fastq", "bar_R1_.fastq"] === (["unp.fastq"], ["bar_R1_.fastq", "bar_R2_.fastq"], ["foo_R1_.fastq", "foo_R2_.fastq"])
groupFastqs fqs = (unpaired, fwd, rev) 
  where
   fwd = matching "_R1_" fqs
   rev = matching "_R2_" fqs 
   unpaired = [x | x <- fqs, not (x `elem` fwd ++ rev)]
   matching str strings = (filter (=~ str) strings) :: [String] 

-- FilePath is an alias for string
-- run pFunc on the paired fastqs, unpFunc on the unpaired fastqs.
-- pFunc runs with forward and reverse.
--mapPairedUnpaired :: (FilePath -> FilePath -> Action ()) -> (FilePath -> Action ()) -> [FilePath] -> Action ()
mapPairedUnpaired pFunc unpFunc fqs = do
  for unp unpFunc
  zipWithM_ pFunc fwd rev  -- like for but zipped
  where (unp, fwd, rev) =  groupFastqs fqs

-- Get all current fastq files, and all sff files which are expected to convert to fastq files.
allFqs = do
 fqs <- getDirectoryFiles "" ["data/*.fastq"]
 sffs <- getDirectoryFiles "" ["data/*.sff"]
 let sffs2fastqs = map (-<.> "fastq") sffs
 return (sffs2fastqs ++ fqs)
  
-- given a string of python code, create a temporary file, write that string to the file,
-- then execute it with the "python" interpreter.
runPython str = do 
  withTempFile $ \file -> do 
    liftIO $ writeFile file fixedStr
    cmd "python" file 
  where
   -- just mangle the text so that the indentation is correct, this allows more readable strings.
   isWhitespace = null . dropWhile (\x -> (x==' ') || (x == '\n'))
   full = dropWhile isWhitespace $ lines str 
   indent = length $ takeWhile (== ' ') $ head full
   fixedStr = unlines $ map (drop indent) full



-- note that shake cannot tell if a rule has changed. so altering want/need and re-running won't work.
       -- () <-  is needed if cmd is not the last statement
