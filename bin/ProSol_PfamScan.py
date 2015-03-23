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



class ProSol_PfamScan:
    
    pfamScanDir     = ScriptsDir.pfamScanDir
    pfam_scan        = pfamScanDir + "pfam_scan.pl "
    
              
    def pfamScan(self, fasta, pfamDir, output ):
        
        try:
            curDirectory = os.getcwd()
            
            os.chdir(ScriptsDir.pfamScanDir)

            call = self.pfam_scan + " -fasta " + fasta + " -dir " + pfamDir + " > " + output
            
            NGS_Util.executeCall(call)
            
            os.chdir(curDirectory)
            
        except Exception:
            
            print traceback.print_exc()
    

