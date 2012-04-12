#! /usr/bin/python
"""
This is a simple script that will take a fasta file,
send it to ncbi for blast over the nr database,
and then store the results in an sqlite database
"""
import argparse
from argparse import RawTextHelpFormatter
import sqlite3
import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def main():
    args = getArgs()
    
    con = setupDb(args.outputdb)
    query_records = getFastaRecords(args.fastafile)
    
    for query_record in query_records:
        blast_record = runBlast(query_record)
        data = getBlastStats(blast_record)
        
        data['query_length'] = len(query_record.upper())
        data['query_id'] = query_record.id

        insertData(con, data)
        
    con.close()
    exit(0)

def getArgs():
    """Use argparse to get command line arguments. Returns"""
    
    parser = argparse.ArgumentParser(description="""
Run BLAST on a givin fasta file and output a summary of the results to an sqlite
database.

The fasta file may contain multiple sequences. For each sequence, blast will be 
run and a row will be added to the database.

Data is stored in the BlastSummary table, with this schema:
    id:             (Integer) Primary Key
    query_id:       (Text)    Query sequence identifier
    query_length:   (Integer) Length of query sequence
    num_alignments: (Integer) Total number of alignments returned
    total_hsps:     (Integer) Total number of hsps across all alignments
    num_hsps:       (Integer) Number of hsps in best alignment
    best_score:     (Real)    Score of best hsp
    best_expect:    (Real)    Expect of best hsp
    best_align_id:  (Text)    Hit identifier of the top scoring alignment
    best_align_def: (Text)    Hit definition of the top scoring alignment
    
The best alignment and best hsp are found by looking for the hsp with the 
highest score and lowest expect, in that order.
"""
    , formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('fastafile', help="Fasta file for query")
    parser.add_argument('outputdb', help="Location of sqlite database.")
    args = parser.parse_args()
    return args

def getFastaRecords(fastafile):
    """Read query file"""
    try:
        records = SeqIO.parse(open(fastafile),  format="fasta")
    except StandardError, e:
        print "There was a problem reading your query file:", e
        exit(1)
    return records

def runBlast(queryRecord):
    """Send fasta query to NCBI and return results"""
    result_handle = NCBIWWW.qblast("blastn",  "nr",  queryRecord.format("fasta"))   
    blast_record = NCBIXML.read(result_handle)
    return blast_record
    
def getBlastStats(blast_record):
    """Gather various stats on the results returned by NCBI"""
    best_align = None
    best_hsp = None
    total_hsp = 0
    
    for alignment in blast_record.alignments:
        if best_align is None:
            best_align = alignment
        for hsp in alignment.hsps:
            total_hsp += 1
            if best_hsp is None:
                best_hsp = hsp
            elif hsp.score > best_hsp.score and hsp.expect < best_hsp.expect:
                best_hsp = hsp
                
    data = {}
    data['total_hsps']     = total_hsp
    data['num_hsps']       = len(best_align.hsps)
    data['num_alignments'] = len(blast_record.alignments)
    data['best_score']     = best_hsp.score
    data['best_expect']    = best_hsp.expect
    data['best_align_id']  = best_align.hit_id
    data['best_align_def'] = best_align.hit_def
    
    return data

def insertData(con, data):
    """Insert a new row to the database"""
    try:
        cur = con.cursor()
        cur.execute("""
INSERT INTO BlastSummary
       (query_id, total_hsps, num_hsps, num_alignments, best_score, best_expect, query_length
       ,best_align_id, best_align_def)
VALUES (:query_id,:total_hsps,:num_hsps,:num_alignments,:best_score,:best_expect,:query_length
       ,:best_align_id,:best_align_def)
""", data)

        con.commit()
    except sqlite3.Error,  e:
        print "Sqlite Error:", e
        sys.exit(1)

def setupDb(outputdb):
    """Create the BlastSummary table if necessary"""
    con = None
    try:
        con = sqlite3.connect(outputdb)
        cur = con.cursor()
        
        cur.execute("""
CREATE TABLE IF NOT EXISTS BlastSummary(
    id             INTEGER PRIMARY KEY AUTOINCREMENT
   ,total_hsps     INTEGER
   ,num_hsps       INTEGER
   ,num_alignments INTEGER
   ,best_score     REAL
   ,best_expect    REAL
   ,query_length   INTEGER
   ,query_id       TEXT
   ,best_align_id  TEXT
   ,best_align_def TEXT)
""")
        con.commit()
    except sqlite3.Error,  e:
        print "Sqlite Error:", e
        sys.exit(1)
        
    return con
    
main()
