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

class ProSol_Hmmer:


    def hmmfetch(self, pfamFile, pfam, output ):
        
        try:
            
            call = "hmmfetch " + pfamFile + " " + pfam + " > " + output
            
            NGS_Util.executeCall(call)
            
        except Exception:
            
            print traceback.print_exc()


    def hmmpress(self, pfamFile):
        
        try:
            
            call = "hmmpress " + pfamFile
            
            NGS_Util.executeCall(call)
            
        except Exception:
            
            print traceback.print_exc()
            
            