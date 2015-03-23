#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""



import os
import sys
import traceback
import NGS_Util


class ProSol_PreProcess:

    orgListFile           = ""  #Name of file containing Organisms List                  #need to be initialized
    orgFastaDir           = ""  #Directory path containing organisms fasta sequences     #need to be initialized         
    accession2speciesFile = ""  #Name of file containing Organisms List with GeneIdentifiers List #need to be initialized
    species2accessionFile = ""  #Name of file containing Organisms List with GeneIdentifiers List #need to be initialized
   

    def initialize(self, orgListFile, orgFastaDir, accession2speciesFile, species2accessionFile ):
    
        try:
            self.orgListFile           = orgListFile
            self.orgFastaDir           = orgFastaDir   
            self.accession2speciesFile = accession2speciesFile
            self.species2accessionFile = species2accessionFile

        except Exception:
            print traceback.print_exc()
    
       
    def create_new_seq_org_list(self):
    
        try:

            print "create_new_seq_org_list"
        
            orgListFile_fh = open(self.orgListFile)
	    accession2speciesFile_fh = open(self.accession2speciesFile,"w")      #output file
	    species2accessionFile_fh = open(self.species2accessionFile,"w")      #output file

            for orgLine in orgListFile_fh:
                
                organismName = orgLine.strip()
                
		print "create_new_seq_org_list: " + organismName

		org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName + ".faa")
		
		org_fasta_fh = open(org_fasta)

		
		for line in org_fasta_fh:

		    if line.startswith(">"):

			seq_id = line.split(">")[1]
			
			if "|" in seq_id:
			    
			    if (">tr|" in line) or (">sp|" in line):

				id = seq_id.split("|")[1].strip()

			    else:
				if " | " in seq_id:
				    id = seq_id.split("|")[0].strip()
				else:
				    id = seq_id.split(" ")[0].strip()


			else:
			    id = seq_id.strip()


			if " " in id:
			    id = id.split(" ")[0].strip()

				    
			accession2speciesFile_fh.write( id + "\t" +  organismName + "\n" )
			species2accessionFile_fh.write( organismName + "\t" +  id + "\n" )
			    
		org_fasta_fh.close()

	    accession2speciesFile_fh.close()
	    species2accessionFile_fh.close()

            orgListFile_fh.close() 
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""

    

    def preProcess(self):
    
        try:

	    self.create_new_seq_org_list()


        except Exception:
            print traceback.print_exc()
        
