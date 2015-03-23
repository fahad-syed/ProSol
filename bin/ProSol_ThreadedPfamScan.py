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
import threading
import NGS_Util
import ScriptsDir
import ProSol_PfamScan

class ProSol_ThreadedPfamScan(threading.Thread):
    
    fasta       = ""
    pfamDir     = ""
    output      = ""
    
    proSol_PfamScan = ProSol_PfamScan.ProSol_PfamScan()

    def __init__(self, fasta, pfamDir, output):
        
        threading.Thread.__init__(self)
        
        self.fasta    = fasta
        self.pfamDir  = pfamDir
        self.output   = output

    def run(self):
        
        self.proSol_PfamScan.pfamScan(self.fasta, self.pfamDir ,self.output)


