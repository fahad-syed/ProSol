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

class ProSol_MCLClustering:


    def makeClusterInputFiles(self, pfamID, abcFile, mclClusteringDir):
        
        try:
            
            mciFile = NGS_Util.createFilePath(mclClusteringDir, pfamID + ".mci")
            tabFile = NGS_Util.createFilePath(mclClusteringDir, pfamID + ".tab")

            call = "mcxload -abc " + abcFile + " --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o " + mciFile + "  -write-tab " + tabFile

            NGS_Util.executeCall(call)
            
            return mciFile, tabFile
            
        except Exception:
            
            print traceback.print_exc()



    def doMCLClustering(self, pfamID, abcFile, mclClusteringDir):
        
        try:

            mciFile, tabFile = self.makeClusterInputFiles(pfamID, abcFile, mclClusteringDir)
            
            I = 1.2
            
            
            for index in range(1,10):
                
                output = NGS_Util.createFilePath(mclClusteringDir, pfamID + ".mci." + str(I).replace(".","") )
                
                call = "mcl " + mciFile + " -I " + str(I) + " -use-tab " + tabFile + " -o " + output

                I += 0.4

                print call
                
                NGS_Util.executeCall(call)
            
        except Exception:
            
            print traceback.print_exc()

