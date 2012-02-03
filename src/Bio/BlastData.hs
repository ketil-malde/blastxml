{- |
   This module implements a hierarchical data structure for BLAST results.

   BLAST is a tool for searching in (biological) sequences for
   similarity.  This library is tested against NCBI-blast version
   2.2.14.  There exist several independent versions of BLAST, so expect some
   incompatbilities if you're using a different BLAST version.

   For parsing BLAST results, the XML format (blastall -m 7) is by far the most
   robust choice, and is implemented in the "Bio.Alignment.BlastXML" module.

   The format is straightforward (and non-recursive).
   For more information on BLAST, check <http://www.ncbi.nlm.nih.gov/Education/BLASTinfo/information3.html>

-}

module Bio.BlastData where

import Data.ByteString.Lazy.Char8 (ByteString)

-- ------------------------------------------------------------

-- | The sequence id, i.e. the first word of the header field.
type SeqId = ByteString

-- | The 'Strand' indicates the direction of the match, i.e. the plain sequence or
--   its reverse complement.
data Strand = Plus | Minus deriving (Read,Show,Eq)

-- | The Aux field in the BLAST output includes match information that depends
--   on the BLAST flavor (blastn, blastx, or blastp).  This data structure captures
--   those variations.
data Aux = Strands !Strand !Strand   -- ^ blastn
         | Frame !Strand !Int      -- ^ blastx
           deriving (Show,Eq)

-- | A 'BlastResult' is the root of the hierarchy.
data BlastResult = BlastResult 
    { blastprogram, blastversion, blastdate :: !ByteString
    , blastreferences :: !ByteString
    , database :: !ByteString
    , dbsequences, dbchars :: !Integer
    , results :: [BlastRecord] }
                   deriving Show

-- | Each query sequence generates a 'BlastRecord'
data BlastRecord = BlastRecord { query :: !SeqId, qlength :: !Int
                               , hits :: [BlastHit] } deriving Show

-- | Each match between a query and a target sequence (or subject)
--   is a 'BlastHit'.
data BlastHit = BlastHit { subject :: !SeqId, slength :: !Int 
                         , matches :: [BlastMatch] } deriving Show
-- | A 'BlastHit' may contain multiple separate matches (typcially when
--   an indel causes a frameshift that blastx is unable to bridge).
data BlastMatch = BlastMatch { bits :: !Double, e_val :: !Double
                             , identity :: (Int,Int)
                             , q_from, q_to, h_from, h_to :: !Int
                             , qseq, hseq :: !ByteString
                             , aux :: !Aux } deriving Show

