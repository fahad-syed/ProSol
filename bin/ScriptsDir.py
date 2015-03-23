#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""


import os
import sys
import traceback
import NGS_Util

projectDir         = "/etc/ProSol/" # needs to be initialized
projectBinDir      = NGS_Util.createDirectoryPath(projectDir,"bin")


##########################################################################################     Blast Toolkit      ##########################################################################################

BlastDir     = "/etc/Blast/ncbi-blast-2.2.28+/bin/" #Needs to be set the user#

BlastDBDir   = NGS_Util.createDirectoryPath(projectDir, "data/BlastDB") 
BlastDustDir = NGS_Util.createDirectoryPath(projectDir, "data/BlastDB")  

############################################################################################################################################################################################################


pfamScanDir = "/etc/PfamScan/"
  

############################################################################################################################################################################################################




ProSol_Thresholding = NGS_Util.createFilePath(projectBinDir,"GetThreshold.py")
