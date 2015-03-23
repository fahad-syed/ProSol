#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""


import os
import sys
import traceback
import NGS_Util
import NGS_Blast

sys.path.append("..")
import ScriptsDir


sys.path.append(ScriptsDir.BLASTScripts)

from buildBlastResult  import combineBlasts


class CreateBlastDB:
   
    orgListFile       = ""  #Name of file containing Organisms List                  #need to be initialized     
    orgFastaDir       = ""  #Directory path containing organisms fasta sequences     #need to be initialized     

    orgBlastDBDir     = ""  #Directory path for organisms BLASTable databases        #need to be initialized     
    orgBlastDustDir   = ""  #Directory path for organisms DUST  files                #need to be initialized     

    ngsBlast          = NGS_Blast.NGS_Blast()

    def initialize(self, orgListFile, orgFastaDir , orgBlastDBDir, orgBlastDustDir):
    
        try:
                   
            self.orgListFile       = orgListFile
            self.orgFastaDir       = orgFastaDir
        
            self.orgBlastDBDir     = orgBlastDBDir
            self.orgBlastDustDir   = orgBlastDustDir
        
        except Exception:
            print traceback.print_exc()
    


    def makeBlastDB(self, organismName):
    
        try:

            print "Make Blast Database: " + organismName
            
            org_fasta    = NGS_Util.createFilePath(self.orgFastaDir, organismName+".faa")
            
            org_dust     = NGS_Util.createFilePath(self.orgBlastDustDir, organismName+"_dust.asnb")

            org_blast_db = NGS_Util.createFilePath(self.orgBlastDBDir, organismName)
            
            if os.path.exists(org_fasta):

                if not os.path.exists(org_blast_db + ".phd") and not os.path.exists(org_blast_db + ".psq"):
                
                    self.ngsBlast.makeProteinBlastDBFromDustFile(org_fasta,org_dust,org_blast_db)
                
                return org_blast_db
            
        except Exception:
            
            print traceback.print_exc()
        
        return ""



    def createOrganismsDB(self):
    
        try:
        
            orgListFile_fh = open(self.orgListFile)

            for line in orgListFile_fh:
                
                organismNameID, organismName = line.strip().split()
                
                self.makeBlastDB(organismName)

            orgListFile_fh.close() 
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""


