#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""

import os
import sys
import traceback
import datetime
import ScriptsDir
import NGS_Util
import ProSol_PreProcess
import ProSol




projectDir         = ScriptsDir.projectDir
projectBinDir      = NGS_Util.createDirectoryPath(projectDir,"bin")
projectDataDir     = NGS_Util.createDirectoryPath(projectDir,"data")
projectResultsDir  = NGS_Util.createDirectoryPath(projectDir,"results")


orgListFile           = NGS_Util.createFilePath(projectDataDir, "org_list.txt")
orgFastaDir           = NGS_Util.createDirectoryPath(projectDataDir, "Genomes_Oct2013")
accession2speciesFile = NGS_Util.createFilePath(projectResultsDir, "accession2species.txt")
species2accessionFile = NGS_Util.createFilePath(projectResultsDir, "species2accession.txt")


interpro_xml          = NGS_Util.createFilePath(projectDataDir, "interpro.xml")
PfamA_hmm_file        = NGS_Util.createFilePath(projectDataDir, "Pfam-A.hmm")

pfamDatfile           = NGS_Util.createFilePath(projectDataDir, "Pfam-A.hmm.dat")
pfamDataDir           = projectDataDir

blastEValue           = 1e-5


def generatePics(currentRunDir):

    try:

        print "generatePics"

        call  =  "cat " +  projectBinDir + "r_analyseClustering.R | R --slave --args " +  currentRunDir
        
        NGS_Util.executeCall(call)
        
        
        picsDir = NGS_Util.createDirectoryPath(currentRunDir,"PicsTables")
        
        epsFile = NGS_Util.createFilePath(picsDir,"SenSpeVsInfAll.eps")
        
        pngFile = NGS_Util.createFilePath(picsDir,"SenSpeVsInfAll.png")
        

        call  =  "convert " +  epsFile + " " +  pngFile
        
        NGS_Util.executeCall(call)

            
    except Exception:
        
        print traceback.print_exc()
        
    return ""



def printHelp():

    try:

        print "\nThe application uses the folloiwng arguments:\n"
        print "     -all                       : use this switch if a user wants to run the application for all of the interopros in the interpro database"
        print "     -i <interpro filename>     : interpro file name"
        print "     -dir  <dir name>           : name of result directory"
        print "     -orgs <org1,org2,...orgn>  : list of organisms seperated by comma. The name of the organisms should be as in <organism name>.faa file name"
        print "     -redoCalc                  : redo calculations for already generated results. If -redoCalc swicth is used then please provide a full result dir path in -dir"
        print "     -resumeAllIPRs             : resume analysis for all of interpros in the interpro database if broken.  If -resumeAllIPRs swicth is used then please provide a full result dir path in -dir"
        print "\n\nAbove, the -i or -all switch is mandatory. A user needs to provide any of these" 

    except Exception:
        
        print traceback.print_exc()
        
    return ""


try:

    iprFileName            = ""
    runDirName             = ""
    useALLInterpros        = False
    organismsToCompareList = []
    displayHelp            = False
    redoCalculations       = False
    resumeALLInterprosCalc = False
    
    
    for index in range(len(sys.argv)):
        
        if sys.argv[index] == "-all":
            useALLInterpros = True
            
        elif sys.argv[index] == "-i":
            iprFileName=sys.argv[index+1]
        
        elif sys.argv[index] == "-dir":
            runDirName = sys.argv[index+1]
            
        elif  sys.argv[index] == "-orgs":
            organismsToCompareList = sys.argv[index+1].split(",")
            
        elif  sys.argv[index] == "-h":
            printHelp()
            displayHelp = True
            
        elif  sys.argv[index] == "-redoCalc":
            redoCalculations = True

        elif  sys.argv[index] == "-resumeAllIPRs":
            resumeALLInterprosCalc = True
    
            
    if not displayHelp:
        
        print "iprFileName            : " + iprFileName
        print "runDirName             : " + runDirName
        print "useALLInterpros        : " + str(useALLInterpros)
        print "organismsToCompareList : " + str(organismsToCompareList) + "\n"
        
        errorFound = False
        

        if (not useALLInterpros) and (len(iprFileName) == 0):
            
            errorFound = True
            print "Error1: Please provide interpro file or use -all switch to run for all interpros."

        elif (useALLInterpros) and (len(iprFileName) > 0):
            
            errorFound = True
            print "Error2: Please provide interpro file or use -all switch to run for all interpros. -all and -i switch cannot be used altogether."


        elif resumeALLInterprosCalc and redoCalculations:
            
            errorFound = True
            print "Error3: -redoCalc and -resumeAllIPRs switches cannnot be used at the same time. Please use any one of them"

        else:
            
            if useALLInterpros  and (len(iprFileName) < 1):
    
                if resumeALLInterprosCalc or redoCalculations:
                
                    if len(runDirName)<1:
                    
                        errorFound = True
                    
                        print "Error4: Please provide an absolute result directory path"
                        
    
                    else:
                        
                        if not os.path.exists(runDirName):
    
                            errorFound = True
                        
                            print "Error5: Directory : "  + runDirName +  " does not exist"
                

            elif (not useALLInterpros) and (len(iprFileName) > 0):
                
                if not os.path.exists(iprFileName):
    
                    errorFound = True
                    
                    print "Error6: Interpro file : "  + iprFileName +  " does not exist"
                    
                    
                if redoCalculations:
                
                    if len(runDirName)<1:
                    
                        errorFound = True
                    
                        print "Error7: Please provide an absolute result directory path"
                        
    
                    else:
                        
                        if not os.path.exists(runDirName):
    
                            errorFound = True
                        
                            print "Error8: Directory : "  + runDirName +  " does not exist"
                            
                            
                if resumeALLInterprosCalc:

                    errorFound = True
                   
                    print "Error9: The -resumeAllIPRs switch is for the whole interpro database"
                            
                              
        if errorFound:
            
            printHelp()
        
        else:
            
            preProcess = ProSol_PreProcess.ProSol_PreProcess()
            preProcess.initialize(orgListFile, orgFastaDir, accession2speciesFile, species2accessionFile)
            preProcess.preProcess()
            
            proSol = ProSol.ProSol()
            
            if not redoCalculations and not resumeALLInterprosCalc:
                            
                strDate = str(datetime.date.isoformat(datetime.date.today()))
        
                if len(runDirName) > 1:
                    currentRunDir = NGS_Util.createDirectoryPath(projectResultsDir, strDate + "_"+ runDirName)                    
                else:                
                    currentRunDir =  NGS_Util.createDirectoryPath(projectResultsDir, strDate)

                    
                NGS_Util.createDirectory(currentRunDir)
                proSol.initialize(orgListFile, orgFastaDir, accession2speciesFile, species2accessionFile, interpro_xml, PfamA_hmm_file, blastEValue, currentRunDir, pfamDatfile, pfamDataDir,organismsToCompareList)

                
                if (useALLInterpros) :
                    
                    proSol.runForAllInterPros()
                    
                else:
                    
                    iprFile_fh = open(iprFileName)
                    
                    proSol.makeOrganismsSequenceIdsDict()
                    
                    for val in iprFile_fh:
                        interproID = val.strip()
                        proSol.runForInterPro(interproID)
                
                    iprFile_fh.close()
                

            elif redoCalculations and not resumeALLInterprosCalc:
                
                if runDirName[len(runDirName)-1] != "/":
                
                    currentRunDir  = runDirName + "/"
                
                else:
                
                    currentRunDir = runDirName
                
                
                
                proSol.initialize(orgListFile, orgFastaDir, accession2speciesFile, species2accessionFile, interpro_xml, PfamA_hmm_file, blastEValue, currentRunDir, pfamDatfile, pfamDataDir,organismsToCompareList)
            
                if (useALLInterpros) :
                    
                    proSol.redoAllInterProsCalculations()
                    
                else:
                    
                    iprFile_fh = open(iprFileName)
                    
                    proSol.makeOrganismsSequenceIdsDict()
                    
                    for val in iprFile_fh:
                        interproID = val.strip()
                        proSol.redoInterproCalculations(interproID)
            
            
            elif not redoCalculations and resumeALLInterprosCalc:

                
                if runDirName[len(runDirName)-1] != "/":
                
                    currentRunDir  = runDirName + "/"
                
                else:
                
                    currentRunDir = runDirName
                
                proSol.initialize(orgListFile, orgFastaDir, accession2speciesFile, species2accessionFile, interpro_xml, PfamA_hmm_file, blastEValue, currentRunDir, pfamDatfile, pfamDataDir,organismsToCompareList)
            
                if (useALLInterpros) :
                    
                    proSol.resumeRunForAllInterPros()

                    

            proSol.calculateNonInformativePfamsInflationValue()
            
            proSol.writeInterproPfamOrganismSelectedOrgsStats()
            
            proSol.writeInterproPfamsMAD()
            
            
            generatePics(currentRunDir)
    
except Exception:
    
    print traceback.print_exc()
