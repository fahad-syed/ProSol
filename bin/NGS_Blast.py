#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""



import os
import sys
import subprocess
import operator
import traceback
import time
import NGS_Util
import ScriptsDir



class NGS_Blast:
    
    blastDir        = ScriptsDir.BlastDir
    numThreads      = " -num_threads 24 "
    segmasker       = blastDir + "segmasker "
    makeblastdb     = blastDir + "makeblastdb "
    blastp          = blastDir + "blastp "
    blast_formatter = blastDir + "blast_formatter "
    

    def makeProteinBlastDBFromDustFile(self, inputFasta, dustFile, blastDB):
        
        try:

            call = self.segmasker + " -in "   + inputFasta + " -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out " + dustFile

            NGS_Util.executeCall(call)

            call = self.makeblastdb + " -in " + inputFasta + " -input_type fasta -dbtype prot -parse_seqids -hash_index -mask_data "   + dustFile + " -out " +  blastDB

            NGS_Util.executeCall(call)
                  
        except Exception:
            print traceback.print_exc()
            
    
    def makeProteinBlastDB(self, inputFasta,blastDB):
        
        try:

            call = self.makeblastdb + " -in " + inputFasta + " -input_type fasta -dbtype prot -parse_seqids -hash_index " + " -out " +  blastDB

            NGS_Util.executeCall(call)
                  
        except Exception:
            print traceback.print_exc()
            
        
    def blastP(self, blastDB, queryFastaFile, outfmt, blastResultFile, eValue ):
        
        try:
            
            call = self.blastp + " -db " + blastDB + " -query " + queryFastaFile + " -outfmt " + str(outfmt)  + " -out " + blastResultFile + " -evalue " + str(eValue) + " " +   self.numThreads

            NGS_Util.executeCall(call)
            
        except Exception:
            
            print traceback.print_exc()
    

    def blastFormatter(self, archive, outfmt, formattedFile):
        	    
        try:
            
            call = self.blast_formatter + " -archive " + archive + " -outfmt " + str(outfmt)  + " -out " + formattedFile

            NGS_Util.executeCall(call)
            
        except Exception:
            
            print traceback.print_exc()


    def blastP_with_Database_Size(self, blastDB, queryFastaFile, outfmt, blastResultFile, eValue, uniprotDBSize ):
        
        try:
            
            call = self.blastp + " -db " + blastDB + " -query " + queryFastaFile + " -outfmt " + str(outfmt)  + " -out " + blastResultFile + " -evalue " + str(eValue) + " " +   self.numThreads + " -dbsize " + str(uniprotDBSize)

            NGS_Util.executeCall(call)
            
        except Exception:
            
            print traceback.print_exc()
    
