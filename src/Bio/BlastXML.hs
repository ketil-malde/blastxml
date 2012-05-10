{-# Language OverloadedStrings #-}

{- | Parse blast XML output.

   If you use a recent version of NCBI BLAST and specify XML output (blastall -m 7),
   this module should be able to parse the result into a hierarchical 'BlastResult'
   structure.

   While the process may consume a bit of memory, the parsing is lazy,
   and file sizes of several gigabytes can be parsed (see e.g. the
   xml2x tool for an example).  To parse XML, we use
   'Text.HTML.TagSoup'. 
-}

module Bio.BlastXML ( readXML
                    , module Bio.BlastData) where

import Bio.BlastData

import qualified Data.ByteString.Lazy.Char8 as B
import Text.HTML.TagSoup
import Control.Monad
import Control.Parallel

type STag = Tag B.ByteString

-- | Parse BLAST results in XML format
readXML :: FilePath -> IO BlastResult
readXML fp = do 
    fc <- B.readFile fp
    when (not $ B.isPrefixOf "<?xml" fc) 
             $ error ("Bio.Sequence.BlastXML.readXML:\n   The file '"
                      ++fp++"' does not look like an XML file - aborting!")
    let ts = parseTags fc
        (h:iters) = breaks (\t -> isTagOpenName "Iteration" t || isTagOpenName "Hit" t) ts
    return $ xml2br h iters

-- | breaks p = groupBy (const (not.p))
breaks :: (a -> Bool) -> [a] -> [[a]]
breaks p (x:xs) = let first = x : takeWhile (not.p) xs
                      rest  = dropWhile (not.p) xs
                  in  rest `par` first : if null rest then [] else breaks p rest
breaks _ []     = []

getFrom :: [STag] -> B.ByteString -> B.ByteString
getFrom list tag = let xs = sections (isTagOpenName tag) list 
                   in if null xs || null (head xs) || (null . drop 1 . head) xs 
                      then error ("Couldn't find tag '"++show tag++"' in\n"++showSome list)
                      else case xs !! 0 !! 1 of 
                             TagText s -> s
                             x -> error ("Unexpeced tag: "++ show x)

-- Use pattern match since 'length' is strict, defeating the purpose.
showSome :: [STag] -> String
showSome a@(_:_:_:_:_:_:_) = (init . show . take 5 $ a)++" ... ]"
showSome a                 = show a

xml2br :: [STag] -> [[STag]] -> BlastResult
xml2br h is = BlastResult { blastprogram = get "BlastOutput_program"
                          , blastversion = bv
                          , blastdate = bd 
                          , blastreferences = get "BlastOutput_reference"
                          , database = get "BlastOutput_db"
                          , dbsequences = 0
                          , dbchars = 0
                          , results = map iter2rec $ breaks (isTagOpenName "Iteration" . head) is
                          }
    where (bv,bd) = B.break (=='[') $ get "BlastOutput_version"
          get = getFrom h

iter2rec :: [[STag]] -> BlastRecord
iter2rec (i:hs) = BlastRecord 
              { query = get "Iteration_query-def"
              , qlength = readI $ get "Iteration_query-len"
              , hits = map hit2hit hs
              }
    where get = getFrom i

iter2rec [] = error "iter2rec: got empty list of sections!"

hit2hit :: [STag] -> BlastHit
hit2hit hs = BlastHit 
             { subject = get "Hit_def"
             , slength = readI $ get "Hit_len"
             , matches = map hsp2match $ partitions (isTagOpenName "Hsp") hs
             }
    where get = getFrom hs


readI :: B.ByteString -> Int
readI x = case B.readInt x of 
  Just (n,_) -> n
  _ -> error ("Couldn't read an Int from string: '"++B.unpack x++"'")

readF :: B.ByteString -> Double
readF = read . B.unpack

hsp2match :: [STag] -> BlastMatch
hsp2match ms = BlastMatch
               { bits   = readF $ get "Hsp_bit-score"
               , e_val  = readF $ get "Hsp_evalue"
               , q_from = readI $ get "Hsp_query-from"
               , q_to   = readI $ get "Hsp_query-to"
               , h_from = readI $ get "Hsp_hit-from"
               , h_to   = readI $ get "Hsp_hit-to"
               , identity = (readI $ get "Hsp_identity", readI $ get "Hsp_align-len")
               , qseq = get "Hsp_qseq"
               , hseq = get "Hsp_hseq"
               -- blastx has query-frame ±1..3 
               -- tblastn has hit-frame
               -- blastn has both hit and query
               -- tblastx has query-frame = 1, hit-frame ±1..3
               , aux = case sections (isTagOpenName "Hsp_hit-frame") ms of
                         [] -> mkFrame $ get "Hsp_query-frame"
                         [(_o:TagText hf:_c)] -> case sections (isTagOpenName "Hsp_query-frame") ms of 
                                                   [] -> mkFrame hf
                                                   [(__o:TagText qf:__c)] -> mkStrands hf qf
                                                   e -> error ("hsp2match: should be tagopen/text/close:\n"++show e)
                         e -> error ("hsp2match: failed to determine frame:\n"++show e)
               }
    where get = getFrom ms
          mkFrame f = Frame (strand' $ signum $ readI f) (abs $ readI f)
          mkStrands h q = Strands (strand' $ readI h) (strand' $ readI q)
          -- ignore frame also for tblastx hits (it can be reconstructed from location)
          strand' :: Int -> Strand
          strand' s = if s > 0 then Plus else Minus


