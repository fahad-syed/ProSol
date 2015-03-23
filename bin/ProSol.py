#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""



import os
import sys
import math
import traceback
import subprocess
import ScriptsDir
import NGS_Util
import DataStrucures_Util
import ProSol_Hmmer
import ProSol_PfamScan
import ProSol_MCLClustering
import ProSol_ThreadedPfamScan
import NGS_Blast
import PfamBlastStatistics
import shutil
from Bio import SeqIO

class ProSol:

    orgListFile           = ""  #Name of file containing Organisms List                  #need to be initialized
    orgFastaDir           = ""  #Directory path containing organisms fasta sequences     #need to be initialized         
    accession2speciesFile = ""  #Name of file containing Organisms List with GeneIdentifiers List #need to be initialized
    species2accessionFile = ""  #Name of file containing Organisms List with GeneIdentifiers List #need to be initialized
    interpro_xml          = ""
    PfamA_hmm_file        = ""
    mockPfamFile          = ""
    orgPfamResultDir      = ""
    pfamBlastSequenceDir  = ""
    blastResultDir        = ""
    blastEValue           = 1e-5
    resultDir             = ""
    interproResultDir     = ""
    mclClusteringDir      = ""
    SensitivitySpecificityDir = ""
    pfamDatfile           = ""
    pfamDataDir           = ""
    
    pfamBlastStaticticsDir=""
    
    interproPfamsDict     = {}
    interproDescDict      = {}
    
    
    clusterPfamOrganismStatsDir = ""
    
    organismsSequenceIDsDict            = {}
    
    interproInflationValuesSelectedDict = {}
    nonInformativeInterproDict          = {}

    organismsToCompareList              = []

    interproPfamOrganismSelectedOrgsStatsDict = {}
    
    interproPfamMADDict = {}
            
            
            
    def initialize(self, orgListFile, orgFastaDir, accession2speciesFile, species2accessionFile, interpro_xml, PfamA_hmm_file, blastEValue, resultDir, pfamDatfile, pfamDataDir, organismsToCompareList):
    
        try:
            self.orgListFile            = orgListFile
            self.orgFastaDir            = orgFastaDir   
            self.accession2speciesFile  = accession2speciesFile
            self.species2accessionFile  = species2accessionFile
            
            self.interpro_xml           = interpro_xml
            self.PfamA_hmm_file         = PfamA_hmm_file
            self.blastEValue            = blastEValue 
            
            self.resultDir              = resultDir
            
            self.pfamDatfile            = pfamDatfile
            self.pfamDataDir            = pfamDataDir
            
            self.organismsToCompareList = organismsToCompareList
            
            self.getAllInterproPfamsAndDescription()
            
        except Exception:
            print traceback.print_exc()


    def makeInterproResultsDirectoriesPath(self, interproID, resultDir):
    
        try:
            self.interproResultDir           = NGS_Util.createDirectoryPath(resultDir,interproID)
            self.mockPfamDir                 = NGS_Util.createDirectoryPath(self.interproResultDir,"MockPfam")
            self.orgPfamResultDir            = NGS_Util.createDirectoryPath(self.interproResultDir,"OrgPFamScan")
            self.mclClusteringDir            = NGS_Util.createDirectoryPath(self.interproResultDir,"MCL_Clustering")
            self.SensitivitySpecificityDir   = NGS_Util.createDirectoryPath(self.interproResultDir,"SensitivitySpecificity")
            self.clusterPfamOrganismStatsDir = NGS_Util.createDirectoryPath(self.interproResultDir,"ClusterPfamOrganismStats")
            self.pfamBlastStaticticsDir      = NGS_Util.createDirectoryPath(self.interproResultDir,"PfamBlastStatictic")

            self.pfamBlastSequenceDir   = self.interproResultDir

            self.blastResultDir         = self.interproResultDir
            
        except Exception:
            print traceback.print_exc()

    def makeInterproResultsDirectories(self, interproID, resultDir):
    
        try:
            self.interproResultDir           = NGS_Util.createDirectory(NGS_Util.createDirectoryPath(resultDir,interproID))
            self.mockPfamDir                 = NGS_Util.createDirectory(NGS_Util.createDirectoryPath(self.interproResultDir,"MockPfam"))
            self.orgPfamResultDir            = NGS_Util.createDirectory(NGS_Util.createDirectoryPath(self.interproResultDir,"OrgPFamScan")) 
            self.mclClusteringDir            = NGS_Util.createDirectory(NGS_Util.createDirectoryPath(self.interproResultDir,"MCL_Clustering"))
            self.SensitivitySpecificityDir   = NGS_Util.createDirectory(NGS_Util.createDirectoryPath(self.interproResultDir,"SensitivitySpecificity"))
            self.clusterPfamOrganismStatsDir = NGS_Util.createDirectory(NGS_Util.createDirectoryPath(self.interproResultDir,"ClusterPfamOrganismStats"))
            self.pfamBlastStaticticsDir      = NGS_Util.createDirectory(NGS_Util.createDirectoryPath(self.interproResultDir,"PfamBlastStatictic"))

            self.pfamBlastSequenceDir   = self.interproResultDir

            self.blastResultDir         = self.interproResultDir
            
        except Exception:
            print traceback.print_exc()
            

            
    def getAllInterproPfamsAndDescription(self):
    
        try:
    
            print "getAllInterproPfams"
        
            interproFile_fh = open(self.interpro_xml)
            
            line = ""
            pfamID = ""
     
            for line in interproFile_fh:
                
                if "<interpro id=\"" in line:
                                    
                    index = line.find("<interpro id=\"") + len("<interpro id=\"")
                    
                    interproID = line[index: line.find("\"",index)]
    
                    interproIDFound = True
                    membersListFound = False
                    
                    
                elif "</interpro>" in line:
    
                        interproIDFound = False
    
    
                if "<name>" in line and interproIDFound:
                    index = line.find("<name>") + len("<name>")
                    
                    name = line[index: line.find("</name>",index)]
                        
                    self.interproDescDict[interproID] = name                
                    
    
    
                if "<member_list>" in line:
                    
                    membersListFound = True
    
                if "</member_list>" in line:
                    membersListFound = False
        
    
      
                if "<db_xref " in line:
                                        
                    if  interproIDFound and membersListFound and  "db=\"PFAM\""  in line:
                        
                       pfamDetailsList = line.strip().split(" ")
                       
                       for item in pfamDetailsList:
                        
                        if "dbkey" in item:
                            pfamID = item.split("\"")[1]
                            
                            if len(pfamID) > 0:
                                DataStrucures_Util.addUniqueListValueToDictionary(interproID, pfamID,self.interproPfamsDict)
                                break
                                
    
            interproFile_fh.close()
                        
            return self.interproPfamsDict, self.interproDescDict
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""
   
            
       
    def getPfamID(self, interproFile, interproID):
    
        try:

            print "get getPfamID: " + interproID
        
            interproFile_fh = open(interproFile)
            
            line = ""
            pfamID = ""
            interproIDFound = False

            for line in interproFile_fh:
                
                if "<interpro id=\"" + interproID + "\"" in line:
                    
                    interproIDFound = True

                elif "</interpro>" in line:
                    
                    interproIDFound = False


                if "<db_xref " in line:
                                        
                    if interproIDFound and "db=\"PFAM\""  in line:
                        
                       pfamDetailsList = line.strip().split(" ")
                       
                       for item in pfamDetailsList:
                        
                        if "dbkey" in item:
                            pfamID = item.split("\"")[1]
                            break
                            

            interproFile_fh.close()
            
            return pfamID
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""




    def grepPfamID_from_PfamA_hmm(self, pfamID, PfamA_hmm_file):
    
        try:

            print "get grepPfamID_from_PfamA_hmm: " + pfamID
            
            fullPfamID = NGS_Util.grep(pfamID,PfamA_hmm_file)            
            
            return fullPfamID[0].split(" ")[3]
        

        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def createMockPfam(self, pfamID, fullPfamID, pfamFile):
    
        try:

            print "createMockPfam: " + fullPfamID 
                        
            self.mockPfamFile = NGS_Util.createFilePath( self.mockPfamDir, "Pfam-A.hmm")
            mockPfamDatFile   = NGS_Util.createFilePath( self.mockPfamDir, "Pfam-A.hmm.dat")
            
            if ( not os.path.exists(mockPfamDatFile) ):

                os.symlink(self.pfamDatfile,mockPfamDatFile)
            
                proSol_Hmmer = ProSol_Hmmer.ProSol_Hmmer()
                
                proSol_Hmmer.hmmfetch(pfamFile, fullPfamID, self.mockPfamFile)
                
                proSol_Hmmer.hmmpress(self.mockPfamFile)
                
                return self.mockPfamFile

        except Exception:
            
            print traceback.print_exc()
            
        return ""
    

    def doOrgPfamScan(self,pfamID ):
    
        try:

            print "doOrgPfamScan"
            
            proSol_PfamScan = ProSol_PfamScan.ProSol_PfamScan()
        
            orgListFile_fh = open(self.orgListFile)

            for orgLine in orgListFile_fh:
                
                organismName = orgLine.strip()
                
                print "pfam scan for: " + organismName
                
                org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName + ".faa")
                
                org_pfamScan = NGS_Util.createFilePath(self.orgPfamResultDir, organismName + "." + pfamID + ".pfamScan")                
                
                proSol_PfamScan.pfamScan(org_fasta, self.mockPfamDir ,org_pfamScan)

            orgListFile_fh.close()
            
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def doThreadedOrgPfamScan(self,pfamID):
    
        try:

            print "doThreadedOrgPfamScan"
            
            threads = []
       
            orgListFile_fh = open(self.orgListFile)

            for orgLine in orgListFile_fh:
                
                organismName = orgLine.strip()
                
                print "pfam scan for: " + organismName
                
                org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName + ".faa")
                
                org_pfamScan = NGS_Util.createFilePath(self.orgPfamResultDir, organismName + "." + pfamID + ".pfamScan")                
                
                
                proSolThreadedPfamScan = ProSol_ThreadedPfamScan.ProSol_ThreadedPfamScan(org_fasta,self.mockPfamDir ,org_pfamScan )
                
                proSolThreadedPfamScan.start()
                
                threads.append(proSolThreadedPfamScan)

            orgListFile_fh.close()

            for pfamScanThread in threads:
                 pfamScanThread.join()
            
            print "finished: doThreadedOrgPfamScan"
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""
     
################################################################################ makeBlast_PfamScanSequence_File -----###################################################################################

    def retrievePfamScan_OrganismIds(self, org_pfamScan_output):
    
        try:

            print "retrievePfamScan_OrganismIds: " + org_pfamScan_output
            
            org_pfamScan_output_fh = open(org_pfamScan_output)
            
            line = ""
            
            seqIdsList=[]

            for line in org_pfamScan_output_fh:

                line = line.strip()

                if not line.startswith("#"):
                    
                    if len(line) > 0:

                        id = line[0:line.index(" ")]
                        seqIdsList.append(id)
                
            org_pfamScan_output_fh.close()
            
            return seqIdsList
     
        except Exception:
            
            print traceback.print_exc()
            
        return []


    def getPfamSearchResultSequences(self, pfamID, organismName, orgPfamScanSeqIdsList, pfamSearchResultSequenceFile):
    
        try:

            print "getPfamSearchResultSequences: " + organismName

        
            pfamId_sequences=[]

            org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName + ".faa")

            org_fasta_fh = open(org_fasta, "rU")
            
            for record in SeqIO.parse(org_fasta_fh, "fasta") :
                
                if record.id in orgPfamScanSeqIdsList:
                    
                    pfamId_sequences.append(record)
                    
            org_fasta_fh.close()
            
                     
            pfamSearchResultSequenceFile_fh = open(pfamSearchResultSequenceFile, "a")
            
            SeqIO.write(pfamId_sequences, pfamSearchResultSequenceFile_fh, "fasta")
            
            pfamSearchResultSequenceFile_fh.close()
            
            
            return pfamSearchResultSequenceFile
        
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def makeBlast_PfamScanSequence_File(self,pfamID,pfamSearchResultSequenceFile):
    
        try:

            print "makeBlast_PfamScanSequence_File"
            
            orgListFile_fh = open(self.orgListFile)

            for orgLine in orgListFile_fh:
                
                organismName = orgLine.strip()
                
                print "Get PfamScan sequences for: " + organismName

                org_pfamScan_output   = NGS_Util.createFilePath(self.orgPfamResultDir, organismName + "." + pfamID + ".pfamScan")
                
                orgPfamScanSeqIdsList = self.retrievePfamScan_OrganismIds(org_pfamScan_output)
                
                self.getPfamSearchResultSequences(pfamID, organismName, orgPfamScanSeqIdsList, pfamSearchResultSequenceFile)

            orgListFile_fh.close()
            
     
        except Exception:
            
            print traceback.print_exc()

###########################################################################----- makeBlast_PfamScanSequence_File #######################################################################################


    def doPfamScanSequenceBlast(self,pfamID, pfamSearchResultSequenceFile):
    
        try:

            print "doPfamScanSequenceBlast"
            
            ngsBlast = NGS_Blast.NGS_Blast()
            
            
            pfamDustFile = NGS_Util.createFilePath(ScriptsDir.BlastDustDir , pfamID + "_dust.asnb")
            pfamBlastDB  = NGS_Util.createFilePath(ScriptsDir.BlastDBDir , pfamID)

            ngsBlast.makeProteinBlastDBFromDustFile(pfamSearchResultSequenceFile,pfamDustFile, pfamBlastDB)
            
            blastResultFile = NGS_Util.createFilePath(self.blastResultDir, pfamID + ".blast")
            
            ngsBlast.blastP(pfamBlastDB, pfamSearchResultSequenceFile, 6, blastResultFile, self.blastEValue )
            
            return blastResultFile 


        except Exception:
            
            print traceback.print_exc()
            
        return ""


################################################################################# getThreshold -----####################################################################################################

    def prepareThresholdingInputFile(self,pfamID, blastResultFile):
    
        try:

            print "prepareThresholdingInputFile"
            
            
            if os.path.exists(blastResultFile):
                
                blastResultFile_fh       = open(blastResultFile)
                thresholdingInputFile    = NGS_Util.createFilePath(self.interproResultDir, pfamID + ".threshold_input")
                thresholdingInputFile_fh = open(thresholdingInputFile, "w")
    
                line = ""            
                blastOutput = []
    
    
                for line in blastResultFile_fh:
                    
                    blastOutput = line.strip().split("\t")
    
                    thresholdingInputFile_fh.write(blastOutput[0] + " " + blastOutput[1]  + " " +  blastOutput[10] + "\n")
                    
    
                blastResultFile_fh.close()
                
                thresholdingInputFile_fh.close()
                
                
                return thresholdingInputFile
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""
            

    def getThreshold(self, pfamID, blastResultFile):
    
        try:

            print "findThresholdingValue"

            if os.path.exists(blastResultFile):
                
                thresholdingInputFile  = self.prepareThresholdingInputFile(pfamID, blastResultFile)
 
                thresholdingOutputFile = NGS_Util.createFilePath(self.interproResultDir, pfamID + ".threshold")
               
                call  =  "python " + ScriptsDir.ProSol_Thresholding + " " + thresholdingInputFile + " > " + thresholdingOutputFile
                
                NGS_Util.executeCall(call)
                
                
                thresholdingOutputFile_fh = open(thresholdingOutputFile, "r")
    
                line = ""            
                threshold = ""
    
                for line in thresholdingOutputFile_fh:
                    
                    if "Threshold Found:" in line:
                        
                        threshold = line.strip().split(":")[1].strip()
                    
                
                print "Threshold Found:" + threshold
                return threshold
            

        except Exception:
            
            print traceback.print_exc()
            
        return ""

###########################################################################----- getThreshold ##########################################################################################################


    def getEValueCutOff(self, threshold):
    
        try:

            print "getEValueCutOff"
            
            thresh = int(threshold)*-1
            
            eValueCutOff = math.pow(10, thresh)
            
            print "EValue CutOFF: " + str(eValueCutOff)
            
            
            return eValueCutOff

        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def makeMclInputFile(self,pfamID, blastResultFile, eValueCutOff):
    
        try:

            print "prepareThresholdingInputFile"
            
            
            if os.path.exists(blastResultFile):
                
                blastResultFile_fh = open(blastResultFile)
                mclInputFile       = NGS_Util.createFilePath(self.mclClusteringDir, pfamID + ".abc")
                mclInputFile_fh    = open(mclInputFile, "w")
    
                line = ""            
                blastOutput = []    
    
                for line in blastResultFile_fh:
                    
                    blastOutput = line.strip().split("\t")
                    
                    if float(blastOutput[10]) <= eValueCutOff:
    
                        mclInputFile_fh.write(blastOutput[0] + "\t" + blastOutput[1]  + "\t" +  blastOutput[10] + "\n")
                    
    
                blastResultFile_fh.close()
                
                mclInputFile_fh.close()
                
                
                return mclInputFile
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def doMCLClustering(self, pfamID, mclInputFile):
    
        try:

            print "doMCLClustering"
            
            mclClustering = ProSol_MCLClustering.ProSol_MCLClustering()
            
            mclClustering.doMCLClustering(pfamID, mclInputFile, self.mclClusteringDir)


        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def doPfamScan(self,pfamID, pfamSearchResultSequenceFile):
    
        try:

            print "doPfamScan"
            
            pfamScan = NGS_Util.createFilePath(self.interproResultDir, pfamID + ".pfamScan")
            
            proSol_PfamScan = ProSol_PfamScan.ProSol_PfamScan()
            
            proSol_PfamScan.pfamScan(pfamSearchResultSequenceFile, self.pfamDataDir ,pfamScan)
            
            return pfamScan

     
        except Exception:
            
            print traceback.print_exc()
            
        return ""
     

################################################################################# Sensitivty Specificity -----########################################################################################

    def makePfamsProteinsDict(self,pfamScan):
        
        try:
            
            print "makePfamsProteinsDict " + pfamScan 

            line = ""
            
            pfamIdsDict={}

            org_pfamScan_output_fh = open(pfamScan)

            for line in org_pfamScan_output_fh:
                
                line =line.strip()
                
                if not line.startswith("#"):
                    
                    if len(line) > 0:
                        
                        pfscanValues = line.split()
                        
                        pfamID = pfscanValues[5] #.split(".")[0]
                        
                        proteinID = pfscanValues[0]
                        
                        pfamIdsDict = DataStrucures_Util.addUniqueListValueToDictionary(pfamID,proteinID,pfamIdsDict)

            org_pfamScan_output_fh.close()
            
            pfamIdsDictFile = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "pfam_protiens_list.txt")
            
            DataStrucures_Util.writeDictionaryListValuesToFile(pfamIdsDict,pfamIdsDictFile)

            return pfamIdsDict
          
        except Exception:
            
            print traceback.print_exc()
            
        return ""
     

    def makeClusterProteinsDict(self, clusterFile):
        
        try:
            
            print "makeClusterProteinsDict " + clusterFile

            line = ""

            counter = 0

            clusterIdsDict={}

            clustert_fh = open(clusterFile)

            for line in clustert_fh:
                
                if len(line) > 0:
                    
                    counter += 1
                    
                    clusterID = "Cluster_" + str(counter)
                    
                    clusterValues = line.strip().split()
                    
                    clusterIdsDict[clusterID] = clusterValues

            clustert_fh.close()
            
            return clusterIdsDict
          
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def calculatePfamClusterSensitivity(self, pfamIdsDict, clusterIdsDict):
        
        try:
            print "calculatePfamClusterSpecificity"
            
            pfamSpecficity = {}
        
            for pfamKey, pfamMembers in pfamIdsDict.iteritems():
                
                pfamClusterSpecificity = []
                
                for clusterKey, clusterMembers in clusterIdsDict.iteritems():
                    
                    pfamMembersCount       = len(pfamMembers)
                    pfamClusterMemberCount = 0
                                        
                    for item in pfamMembers:
                        
                        if clusterMembers.count(item) > 0:
                            
                            pfamClusterMemberCount +=1 

                    specificity= (pfamClusterMemberCount* 1.0)/pfamMembersCount
                    pfamClusterSpecificity.append(specificity)
                    
                
                pfamSpecficity[pfamKey] = pfamClusterSpecificity
                
            return pfamSpecficity
                
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def calculatePfamClusterSpecificity(self, pfamIdsDict, clusterIdsDict):
        
        try:
            
            print "calculatePfamClusterSensitivity"

            pfamSensitivity = {}
        
            for pfamKey, pfamMembers in pfamIdsDict.iteritems():

                pfamClusterSensitivity = []
                
                for clusterKey, clusterMembers in clusterIdsDict.iteritems():
                    
                    clusterMembersCount    = len(clusterMembers)
                    pfamClusterMemberCount = 0
                    
                    for item in pfamMembers:
                        
                        if clusterMembers.count(item) > 0:
                            
                            pfamClusterMemberCount +=1 
                    
                    sensitivity = (pfamClusterMemberCount*1.0)/clusterMembersCount
                    pfamClusterSensitivity.append(sensitivity)
                            
                
                pfamSensitivity[pfamKey] = pfamClusterSensitivity
                
            return pfamSensitivity
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def calculatePfamClusterPercentAverageValues(self, pfamClusterValues):
        
        try:
            
            print "calculatePfamClusterPercentAverageValues"
            
            percentAverage = {}
        
            for pfamKey, pfamClusterValuesList in pfamClusterValues.iteritems():
                
                sumClusterValues = 0.0
                
                for clusterValue in pfamClusterValuesList:
                    
                    sumClusterValues += clusterValue
                    
                    
                percentAverage[pfamKey] = round( (sumClusterValues/len(pfamClusterValues))*100 , 2 )
                
            return percentAverage
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def calculatePfamClusterSensitivityPercentAverageValues(self, pfamSensitivity):
        
        try:
            
            print "calculatePfamClusterSensitivityPercentAverageValues"
            
            percentAverage = {}
        
            for pfamKey, pfamSensitivityValues in pfamSensitivity.iteritems():
                
                maxSensitivityValue = max(pfamSensitivityValues) * 100.00
                    
                percentAverage[pfamKey] = maxSensitivityValue
                
            return percentAverage
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def calculatePfamClusterSpecificityPercentAverageValues(self, pfamSpecficity, pfamSensitivity):
        
        try:
            
            print "calculatePfamClusterSensitivityPercentAverageValues"
            
            percentAverage = {}
        
            for pfamKey, pfamSensitivityValues in pfamSensitivity.iteritems():
                
                sumValues = 0.0
                
                pfamSpecificityValues = pfamSpecficity[pfamKey]
                
                pfamSpecificityValues_by_pfamSensitivityValues = map(lambda x,y: x*y, pfamSpecificityValues,pfamSensitivityValues)

                
                for value in pfamSpecificityValues_by_pfamSensitivityValues:
                    
                    sumValues += value
                    
                    
                percentAverage[pfamKey] = sumValues*100.00
                
            return percentAverage
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def calculateAllPfamsSpecificitySensitivity(self, pfamScan, pfamID):
        
        try:
            print "calculateAllPfamsSpecificitySensitivity"
            
            inflationClustersDict = {}
            
            pfamIdsDict = self.makePfamsProteinsDict(pfamScan)
            
            
            print "--------------pfamIdsDict"
            print pfamIdsDict
            
            
            self.writePfamStatistics(pfamIdsDict)

            clusterStatsFile = NGS_Util.createFilePath(self.mclClusteringDir, "cluster_stats.txt")
            clusterStatsFile_fh = open(clusterStatsFile, "w")
            clusterStatsFile_fh.write("InflationValue\tTotalNumberOfClusters\tTotalNumberOfOrphanClusters\n")
            clusterStatsFile_fh.close()


            inflationValue = 1.2
            
            for index in range(1,10):
                
                strInflationValue                       = str(inflationValue).replace(".","") #str(inflationValue).replace(".","")
                
                clusterFile                             = NGS_Util.createFilePath(self.mclClusteringDir, pfamID + ".mci." + strInflationValue )
                
                clusterIdsDict                          = self.makeClusterProteinsDict(clusterFile)
                
                pfamSpecficity                          = self.calculatePfamClusterSpecificity(pfamIdsDict, clusterIdsDict)
                
                pfamSensitivity                         = self.calculatePfamClusterSensitivity(pfamIdsDict, clusterIdsDict)
               
#                pfamPercentAverageSpecficity            = self.calculatePfamClusterPercentAverageValues(pfamSpecficity)
                
#                pfamPercentAverageSensitivity           = self.calculatePfamClusterPercentAverageValues(pfamSensitivity)


                pfamPercentAverageSpecficity            = self.calculatePfamClusterSpecificityPercentAverageValues(pfamSpecficity, pfamSensitivity)
                
                pfamPercentAverageSensitivity           = self.calculatePfamClusterSensitivityPercentAverageValues(pfamSensitivity)

                
                pfamPercentAverageSpecficitySensitivity = DataStrucures_Util.mergeDictionaries(pfamPercentAverageSensitivity,pfamPercentAverageSpecficity)
            
                inflationClustersDict[strInflationValue]= [pfamIdsDict,clusterIdsDict,pfamSpecficity,pfamSensitivity,pfamPercentAverageSpecficitySensitivity]
            
                inflationValue += 0.4
                

                self.writeClusterStatistics(strInflationValue, clusterIdsDict)

                
                pfamClusterSpecificityFile                  = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "pfam_cluster_" + strInflationValue + ".specificity")
                pfamClusterSensitivityFile                  = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "pfam_cluster_" + strInflationValue + ".sensitivity")
                pfamPercentAverageSpecficitySensitivityFile = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "pfam_cluster_" + strInflationValue + ".PercentAverageSpecficitySensitivity")
            
            
                clusterIdsHeader = " \t" + DataStrucures_Util.getDictionaryKeysAsTabSeperatedValues(clusterIdsDict)
                
                DataStrucures_Util.writeDictionaryListValuesToFileWithHeader(clusterIdsHeader,pfamSpecficity,pfamClusterSpecificityFile)                
                
                DataStrucures_Util.writeDictionaryListValuesToFileWithHeader(clusterIdsHeader,pfamSensitivity,pfamClusterSensitivityFile)
                
                DataStrucures_Util.writeDictionaryListValuesToFileWithHeader("PfamID\tSensitivity\tSpecficity", pfamPercentAverageSpecficitySensitivity, pfamPercentAverageSpecficitySensitivityFile)

            return inflationClustersDict
                
        except Exception:
            
            print traceback.print_exc()
            
        return ""

## fix variance prob
    def calculateBestPfamVariance(self, pfamID, fullPfamID, inflationClustersDict):
        
        try:
            
            print "calculatePfamSpecificitySensitivity"
            
            pfamPercentAverage                        = {}
            mostVariedPfamList                        = []
            pfamVarianceDict                          = {}
            pfamInflationVarianceDict                 = {}
            mostVariedPfams_vs_secondMostAbundantPfam = {}
            
            for inflationValue, values in inflationClustersDict.iteritems():

                pfamIdsDict                             = values[0]
                
                clusterIdsDict                          = values[1]
                
                pfamSpecficity                          = values[2]
                
                pfamSensitivity                         = values[3]
                                
                pfamPercentAverageSpecficitySensitivity = values[4]
                
                mostVariedPfamInflationID, pfamClusterHits, pfamNormalizedClusterHits, pfamClusterHitsMean, pfamClusterHitsVariance = self.getMostVariedPfam(clusterIdsDict,pfamIdsDict)
                
                pfamVarianceDict[inflationValue] = [mostVariedPfamInflationID, pfamClusterHits, pfamNormalizedClusterHits, pfamClusterHitsMean, pfamClusterHitsVariance]

                pfamInflationVarianceDict[inflationValue] = [mostVariedPfamInflationID, pfamClusterHitsVariance[mostVariedPfamInflationID]]
               
                pfamClusterHitsStatsFile = NGS_Util.createFilePath(self.SensitivitySpecificityDir, str(inflationValue) + "_Pfam_Cluster_Hits_Variance.txt")
                self.writePfamClusterHitsStatistics(pfamNormalizedClusterHits, pfamClusterHitsMean, pfamClusterHitsVariance,pfamClusterHitsStatsFile)


            pfamsIflationVarianceFile = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "most_varied_pfam_per_inflation_value.txt")
            DataStrucures_Util.writeDictionaryListValuesToFileWithHeader("InflationValue\tPfam\tVariance",pfamInflationVarianceDict,pfamsIflationVarianceFile)                

            
            mostVariedPfamVariance = -1.0
            
            
            for inflationValue in sorted(pfamInflationVarianceDict.keys()):
                
                mostVariedPfamInflationID = pfamInflationVarianceDict[inflationValue][0]               
                pfamClusterHitsVariance   = pfamInflationVarianceDict[inflationValue][1]

                if pfamClusterHitsVariance > mostVariedPfamVariance:

                    mostVariedPfamID          = mostVariedPfamInflationID                    
                    mostVariedPfamVariance    = pfamClusterHitsVariance
                    maxVarianceInflationValue = inflationValue
                    
                    
            mostVariedPfamList = [mostVariedPfamID,mostVariedPfamVariance, maxVarianceInflationValue]
                
                
            secondMostAbundantPfam         = self.getSecondMostAbundantPfam(pfamIdsDict)
            secondMostAbundantPfamCount    = len(pfamIdsDict[secondMostAbundantPfam])
#            secondMostAbundantPfamVariance = pfamInflationVarianceDict[secondMostAbundantPfam]
                
            mostVariedPfamCount            = len(pfamIdsDict[mostVariedPfamID])
            mostVariedPfamVariance         = mostVariedPfamVariance

                
#            mostVariedPfams_vs_secondMostAbundantPfam[inflationValue] = [secondMostAbundantPfam,secondMostAbundantPfamCount,secondMostAbundantPfamVariance,mostVariedPfamID,mostVariedPfamCount,mostVariedPfamVariance]
            mostVariedPfams_vs_secondMostAbundantPfam[inflationValue] = [secondMostAbundantPfam,secondMostAbundantPfamCount, mostVariedPfamID,mostVariedPfamCount,mostVariedPfamVariance]

            for inflationValue, values in inflationClustersDict.iteritems():

                pfamPercentAverageSpecficitySensitivity = values[4]            

                for pfamKey, pfamPercentAverageSpecficitySensitivityValues in pfamPercentAverageSpecficitySensitivity.iteritems():
                    
                    if pfamKey == mostVariedPfamID:
                        
                        pfamPercentAverage[inflationValue] = pfamPercentAverageSpecficitySensitivityValues
            
            
            pfamPercentAverageSpecficitySensitivityFile = NGS_Util.createFilePath(self.SensitivitySpecificityDir, mostVariedPfamID + ".PercentAverageSpecficitySensitivity")
            DataStrucures_Util.writeDictionaryListValuesToFileWithHeader("InflationValue\tSensitivity\tSpecficity",pfamPercentAverage,pfamPercentAverageSpecficitySensitivityFile)                
                
            mostVariedPfamsFile = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "mostVariedPfams.txt")
            DataStrucures_Util.writeListValuesToFileWithHeader("Pfam\tVariance\tInflationValue",mostVariedPfamList,mostVariedPfamsFile)                
                
            mostVariedPfams_vs_secondMostAbundantPfam_File = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "mostVariedPfams_vs_secondMostAbundantPfam.txt")
#            DataStrucures_Util.writeDictionaryListValuesToFileWithHeader("InflationValue\tsecondMostAbundantPfam\tsecondMostAbundantPfamCount\tsecondMostAbundantPfamVariance\tmostVariedPfamID\tmostVariedPfamCount\tmostVariedPfamVariance",mostVariedPfams_vs_secondMostAbundantPfam,mostVariedPfams_vs_secondMostAbundantPfam_File)
            DataStrucures_Util.writeDictionaryListValuesToFileWithHeader("InflationValue\tsecondMostAbundantPfam\tsecondMostAbundantPfamCount\tmostVariedPfamID\tmostVariedPfamCount\tmostVariedPfamVariance",mostVariedPfams_vs_secondMostAbundantPfam,mostVariedPfams_vs_secondMostAbundantPfam_File)
            
            
            return mostVariedPfamList, pfamVarianceDict

            
        except Exception:
            
            print traceback.print_exc()
        
            return ""


    def doSpecificitySensitivityCalculations(self, pfamScan, pfamID, fullPfamID):
        
        try:
            
            print "doSpecificitySensitivityCalculations"
            
            inflationClustersDict = self.calculateAllPfamsSpecificitySensitivity(pfamScan, pfamID)
            
            print "-------------inflationClustersDict"
            print inflationClustersDict
            
            mostVariedPfamList, pfamVarianceDict = self.calculateBestPfamVariance(pfamID, fullPfamID, inflationClustersDict)
            
            
            return inflationClustersDict, mostVariedPfamList, pfamVarianceDict
            
        except Exception:
            
            print traceback.print_exc()


    def doInflationValueCalculations(self, inflationClustersDict, mostVariedPfamList):
        
        try:
            
            print "doSpecificitySensitivityCalculations"
            
            sensivitySpecfivityDiffDict = self.getSensivitySpeficityUnitVector(mostVariedPfamList, inflationClustersDict)
            
            mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff = self.findBestInflationValue(sensivitySpecfivityDiffDict)
            
            return mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff

        except Exception:
            
            print traceback.print_exc()



###########################################################################----- Sensitivty Specificity###############################################################################################



###########################################################################Variance Calculations-------###############################################################################################


    def createPfamClusterDS(self, clusterIdsDict, pfamIdsDict):
        
        try:
            
            print "createPfamClusterDS"
            
            pfamClusterHits = {}
            
            for pfamID in pfamIdsDict.iterkeys():

                if not pfamClusterHits.has_key(pfamID):
                    
                    pfamClusterHits[pfamID] = {}
                
                for clusterID in clusterIdsDict.iterkeys():
                    
                    if not pfamClusterHits[pfamID].has_key(clusterID):
                            
                        pfamClusterHits[pfamID][clusterID] = 0
                
            return pfamClusterHits
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""
           
            
    def getPfamClusterHits(self, clusterIdsDict, pfamIdsDict):
        
        try:
            
            print "getPfamClusterHits"
            
            pfamClusterHits = {}
            
            pfamClusterHits = self.createPfamClusterDS(clusterIdsDict, pfamIdsDict)
            
        
            for pfamID, pfamProteinsList in pfamIdsDict.iteritems():
                
                for protein in pfamProteinsList:
                    
                    for clusterID, clusterIdsList in clusterIdsDict.iteritems():
                        
                        if protein in clusterIdsList:
                            
                            pfamClusterHits[pfamID][clusterID] += 1
                    
            return pfamClusterHits
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""

            
    def getPfamNormalizedClusterHits(self, clusterIdsDict, pfamClusterHits):
        
        try:
            
            print "getPfamNormalizedClusterHits"
            
        
            for pfamID, clusterHits in pfamClusterHits.iteritems():
                
                for clusterID in clusterHits.iterkeys():
                    
                    pfamClusterHits[pfamID][clusterID] = (pfamClusterHits[pfamID][clusterID]*1.0)/len(clusterIdsDict[clusterID])
                    
            return pfamClusterHits
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def getPfamClusterHitsMean(self, pfamClusterHits):
        
        try:
            
            print "getPfamClusterHits"
            
            pfamClusterHitsMean = {}
            
            for pfamID, clusterHits in pfamClusterHits.iteritems():
                
                sumOfCounts = 0
                
                for clusterID in clusterHits.iterkeys():
                    
                    sumOfCounts += clusterHits[clusterID]
                
                pfamClusterHitsMean[pfamID] = (sumOfCounts*1.0)/len(clusterHits)
                    
            return pfamClusterHitsMean
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""
    
    
    def getPfamClusterHitsVariance(self, pfamClusterHits, pfamClusterHitsMean):
        
        try:
            
            print "getPfamClusterHitsVariance"
            
            pfamClusterHitsVariance = {}
            
            for pfamID, clusterHits in pfamClusterHits.iteritems():
                
                sumOfPfamClusterHitsVariance = 0
                
                for clusterID in clusterHits.iterkeys():
                    
                    sumOfPfamClusterHitsVariance += math.pow(clusterHits[clusterID] - pfamClusterHitsMean[pfamID], 2) 
                
                pfamClusterHitsVariance[pfamID] = round( sumOfPfamClusterHitsVariance/len(clusterHits) , 5 )
                    
            return pfamClusterHitsVariance

        except Exception:
            
            print traceback.print_exc()
            
        return ""

    
    def getMostVariedPfam(self, clusterIdsDict, pfamIdsDict):
        
        try:
            
            print "getMostVariedPfam"

            
            pfamClusterHits           = self.getPfamClusterHits(clusterIdsDict, pfamIdsDict)
            
            pfamNormalizedClusterHits = self.getPfamNormalizedClusterHits(clusterIdsDict, pfamClusterHits)
            
            pfamClusterHitsMean       = self.getPfamClusterHitsMean(pfamNormalizedClusterHits) ##self.getPfamClusterHitsMean(pfamClusterHits)
            
            pfamClusterHitsVariance   = self.getPfamClusterHitsVariance(pfamNormalizedClusterHits, pfamClusterHitsMean) ##self.getPfamClusterHitsVariance(pfamClusterHits, pfamClusterHitsMean)
            
            maxVariance = 0.0
            
            mostVariedPfamID = ""
            
            for pfamID in pfamClusterHitsVariance:
                
                if pfamClusterHitsVariance[pfamID] >= maxVariance:
                    
                    maxVariance = pfamClusterHitsVariance[pfamID]
                    mostVariedPfamID = pfamID
            
            return mostVariedPfamID, pfamClusterHits, pfamNormalizedClusterHits, pfamClusterHitsMean, pfamClusterHitsVariance
        
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""

                    
    def getSecondMostAbundantPfam(self, pfamIdsDict):
        
        try:
            print "getSecondMostAbundantPfam"

            
            pfamIdsCountDict={}
            
            for pfamKey, pfamMembers in pfamIdsDict.iteritems():
                
                count = len(pfamMembers)
                
                pfamIdsCountDict[pfamKey]=count

            sorted_pfamIdKeys = sorted(pfamIdsCountDict, key=pfamIdsCountDict.get, reverse=True)
            
            if len(sorted_pfamIdKeys) == 1:
                return sorted_pfamIdKeys[0]
                
            return sorted_pfamIdKeys[1]

            
        except Exception:
            
            print traceback.print_exc()
            

###########################################################################----- Variance Calculations ###############################################################################################



########################################################################### Inflation Value Selection-------##########################################################################################

    def getSensivitySpeficityNorm(self, mostVariedPfamList, inflationClustersDict):
        
        try:
            
            print "getSensivitySpeficityNorm"
            
            specficity = 0
            sensitivity = 0
            
            mostVariedPfamID = mostVariedPfamList[0]
            
            for inflationValue in inflationClustersDict.iterkeys():
                
                inflationCluster = inflationClustersDict[inflationValue]
                
                pfamPercentAverageSpecficitySensitivity = inflationCluster[4]
                
                specficity  += math.pow(pfamPercentAverageSpecficitySensitivity[mostVariedPfamID][0],2)
                sensitivity += math.pow(pfamPercentAverageSpecficitySensitivity[mostVariedPfamID][1],2)
                
            specficityNorm = math.sqrt(specficity)
            sensitivityNorm = math.sqrt(sensitivity)

            
            return specficityNorm, sensitivityNorm

            
        except Exception:
            
            print traceback.print_exc()
        
            return ""


    def getSensivitySpeficityUnitVector(self, mostVariedPfamList, inflationClustersDict):
        
        try:
            
            print "getSensivitySpeficityUnitVector"
            
            sensivitySpecfivityDiffDict={}
            
            specficityNorm, sensitivityNorm =  self.getSensivitySpeficityNorm(mostVariedPfamList, inflationClustersDict)
            

            mostVariedPfamID = mostVariedPfamList[0]
            
            for inflationValue in inflationClustersDict.iterkeys():
    
                inflationCluster = inflationClustersDict[inflationValue]
                
                pfamPercentAverageSpecficitySensitivity = inflationCluster[4]
                
                specficity  = pfamPercentAverageSpecficitySensitivity[mostVariedPfamID][0]/specficityNorm
                sensitivity = pfamPercentAverageSpecficitySensitivity[mostVariedPfamID][1]/sensitivityNorm
                
                specficity_sensitivity_diff = abs(specficity - sensitivity)
                
                sensivitySpecfivityDiffDict[inflationValue] = [mostVariedPfamID,specficity,sensitivity,specficity_sensitivity_diff  ]
                
            best_pfam_specficity_sensitivity_diff_File = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "best_pfam_specficity_sensitivity_diff.txt")
            DataStrucures_Util.writeDictionaryListValuesToFileWithHeader("MostVariedPfamID\tspecficity\tsensitivity\tspecficity sensitivity difference",sensivitySpecfivityDiffDict,best_pfam_specficity_sensitivity_diff_File)                
            

            return sensivitySpecfivityDiffDict

            
        except Exception:
            
            print traceback.print_exc()
        
            return ""


    def selectInflationValue(self, pfamVarianceDict, inflationClustersDict):
        
        try:
            
            print "select Inflation Value"
            
            sensivitySpecfivityDiffDict={}
            
            for inflationValue, values in pfamVarianceDict.iteritems():

                mostVariedPfamID = values[0]
                
                inflationCluster = inflationClustersDict[inflationValue]
                
                pfamPercentAverageSpecficitySensitivity = inflationCluster[4]
                
                specficity  = pfamPercentAverageSpecficitySensitivity[mostVariedPfamID][0]
                sensitivity = pfamPercentAverageSpecficitySensitivity[mostVariedPfamID][1]
                
                specficity_sensitivity_diff = abs(specficity - sensitivity)
                
                sensivitySpecfivityDiffDict[inflationValue] = [mostVariedPfamID,specficity,sensitivity,specficity_sensitivity_diff  ]
                
            best_pfam_specficity_sensitivity_diff_File = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "best_pfam_specficity_sensitivity_diff.txt")
            DataStrucures_Util.writeDictionaryListValuesToFileWithHeader("MostVariedPfamID\tspecficity\tsensitivity\tspecficity sensitivity difference",sensivitySpecfivityDiffDict,best_pfam_specficity_sensitivity_diff_File)                
            
            return sensivitySpecfivityDiffDict

            
        except Exception:
            
            print traceback.print_exc()
        
            return ""


    def findBestInflationValue(self, sensivitySpecfivityDiffDict):
        
        try:
            
            print "findBestInflationValue"
            
            
            inflationValuesList = sorted(sensivitySpecfivityDiffDict)

            bestInflationValue  = inflationValuesList[0]

            mostVariedPfamID    = sensivitySpecfivityDiffDict[bestInflationValue][0]            
            
            min_specficity_sensitivity_diff = sensivitySpecfivityDiffDict[bestInflationValue][3]
            
#            for inflationValue, values in sensivitySpecfivityDiffDict.iteritems():

            for inflationValue  in inflationValuesList:
                
                specficity_sensitivity_diff = sensivitySpecfivityDiffDict[inflationValue][3]
                
                if specficity_sensitivity_diff < min_specficity_sensitivity_diff :
                    
                    mostVariedPfamID   = sensivitySpecfivityDiffDict[inflationValue][0]
                    bestInflationValue = inflationValue
                    min_specficity_sensitivity_diff = specficity_sensitivity_diff
                    
            

            pfamStatsFile = NGS_Util.createFilePath(self.interproResultDir, "best_pfam_inflation_value.txt")
            pfamStatsFile_fh = open(pfamStatsFile, "w")
            
            pfamStatsFile_fh.write( "PfamID\tBestInflationValue\tMin_Specficity_Sensitivity_Difference\n")
            pfamStatsFile_fh.write( mostVariedPfamID + "\t" + str(float(bestInflationValue)/10) + "\t" + str(min_specficity_sensitivity_diff) + "\n")
                
            pfamStatsFile_fh.close()
            
            
            return mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff

            
        except Exception:
            
            print traceback.print_exc()
        
            return ""


###########################################################################----- Inflation Value Selection############################################################################################



########################################################################### Cluster Pfam Organism Statistics-------###################################################################################


    def makeOrganismsSequenceIdsDict(self):
    
        try:

            print "makeOrganismSequenceIdsDict:"

            species2accessionFile_fh = open(self.species2accessionFile)
            
            values = []
            line = ""
            
            for line in species2accessionFile_fh:
                
                if len(line) > 1:

                    values = line.split("\t")
                    
                    self.organismsSequenceIDsDict = DataStrucures_Util.addListValueToDictionary(values[0].strip(), values[1].strip(),self.organismsSequenceIDsDict)
        
        except Exception:
            
            print traceback.print_exc()
            


    def createClusterPfamOrganismDS(self, clusterIdsDict, pfamIdsDict):
        
        try:
            
            print "createPfamClusterDS"
            
            clusterPfamOrganismHits = {}
            
            for clusterID in clusterIdsDict.iterkeys():
                
                if not clusterPfamOrganismHits.has_key(clusterID):
                    
                    pfamHitsPercentDict = {}
                    organismHitsPercentDict  = {}
                    
                    for pfamID in pfamIdsDict.iterkeys():
                    
                        if not pfamHitsPercentDict.has_key(pfamID):
                                
                            pfamHitsPercentDict[pfamID] = 0.0
                            
                    for organismID in self.organismsSequenceIDsDict.iterkeys():
                    
                        if not organismHitsPercentDict.has_key(organismID):
                                
                            organismHitsPercentDict[organismID] = 0.0

                    clusterPfamOrganismHits[clusterID] = [pfamHitsPercentDict, organismHitsPercentDict]
                                        
            return clusterPfamOrganismHits
         
        except Exception:
            
            print traceback.print_exc()
            
        return ""



    def calculateClusterPfamOrganismStatistics(self,inflationClustersDict):
        
        try:
            
            print "calculateClusterPfamOrganismStatistics"

            clusterPfamOrganismHits = {}
                        
            for inflationValue, values in inflationClustersDict.iteritems():

                pfamIdsDict             = values[0]
                clusterIdsDict          = values[1]
                
                clusterPfamOrganismHits = self.createClusterPfamOrganismDS(clusterIdsDict, pfamIdsDict)
                
                for clusterID, clusterHits in clusterIdsDict.iteritems():
                    
                    pfamHitsPercentDict     = clusterPfamOrganismHits[clusterID][0]
                    organismHitsPercentDict = clusterPfamOrganismHits[clusterID][1]
                    
                    for pfamID, pfamProteins in pfamIdsDict.iteritems():
                        
                        for protein in pfamProteins:
                            
                            if protein in clusterHits:
                                
                                pfamHitsPercentDict[pfamID] += 1.0

                    for pfamID, pfamProteins in pfamIdsDict.iteritems():
                        
                        pfamHitsPercentDict[pfamID] = (pfamHitsPercentDict[pfamID] * 100.0) / len(pfamProteins)


                                
                    for organism, organismProteinsList in self.organismsSequenceIDsDict.iteritems():
                                                    
                        for protein in organismProteinsList:

                            if protein in clusterHits:
                            
                                organismHitsPercentDict[organism] += 1.0

                    for organism, organismProteinsList in self.organismsSequenceIDsDict.iteritems():
                                                                               
                            organismHitsPercentDict[organism] = (organismHitsPercentDict[organism] * 100.0) / len(organismProteinsList)

                            
                clusterPfamOrganismHitsFile = NGS_Util.createFilePath(self.clusterPfamOrganismStatsDir, "ClusterPfamOrganismStats_" + inflationValue + ".txt")
                self.writeClusterPfamOrganismStatistics(clusterPfamOrganismHits, clusterPfamOrganismHitsFile)
            
            
            return clusterPfamOrganismHits

            
        except Exception:
            
            print traceback.print_exc()
        
            return ""



###########################################################################----- Cluster Pfam Organism Statistics####################################################################################



########################################################################### Write to File-------#####################################################################################################


    def writePfamStatistics(self, pfamIdsDict):
        
        try:
            print "writePfamStatistics"

            pfamStatsFile = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "pfam_protiens_list.stats")
            pfamStatsFile_fh = open(pfamStatsFile, "w")
            
            
            for pfamKey, pfamMembers in pfamIdsDict.iteritems():
                
                count = len(pfamMembers)
 
                pfamStatsFile_fh.write( pfamKey + "\t" + str(count) + "\n")
                
            pfamStatsFile_fh.close()
            
        except Exception:
            
            print traceback.print_exc()


    def writeClusterStatistics(self, inflationValue, clusterIdsDict):
        
        try:
            
            print "writeClusterStatistics"
            
            totalNumberOfClusters = len(clusterIdsDict)
            
            totalNumberOfOrphanClusters = 0
            
            for clusterKey, clusterMembers in clusterIdsDict.iteritems():
                
                count = len(clusterMembers)
                    
                if count == 1:
                     totalNumberOfOrphanClusters +=1 
 
                
            clusterStatsFile = NGS_Util.createFilePath(self.mclClusteringDir, "cluster_stats.txt")
            clusterStatsFile_fh = open(clusterStatsFile, "a")
            
            clusterStatsFile_fh.write( str(inflationValue) + "\t" + str(totalNumberOfClusters) + "\t" + str(totalNumberOfOrphanClusters) + "\n")
            clusterStatsFile_fh.close()
            
        except Exception:
            
            print traceback.print_exc()


    def writePfamClusterHits(self, pfamIdsDict):
        
        try:
            print "writePfamStatistics"

            pfamStatsFile = NGS_Util.createFilePath(self.SensitivitySpecificityDir, "pfam_protiens_list.stats")
            pfamStatsFile_fh = open(pfamStatsFile, "w")
            
            
            for pfamKey, pfamMembers in pfamIdsDict.iteritems():
                
                count = len(pfamMembers)
 
                pfamStatsFile_fh.write( pfamKey + "\t" + str(count) + "\n")
                
            pfamStatsFile_fh.close()
            
        except Exception:
            
            print traceback.print_exc()
            
            
    def writePfamClusterHitsStatistics(self, pfamClusterHits, pfamClusterHitsMean, pfamClusterHitsVariance, pfamClusterHitsStatsFile ):
        
        try:
            
            print "writePfamClusterHitsStatistics"
            
            pfamClusterHitsStatsFile_fh = open(pfamClusterHitsStatsFile, "w")
            
            header = "Pfam Ids"
            

            clusterHits = pfamClusterHits.items()[0][1]
            
            for clusterID in clusterHits.iterkeys():

                header += "\t" + str(clusterID)

            
            
            header += "\tMean\tVariance\n"
            
            pfamClusterHitsStatsFile_fh.write(header)

            for pfamID, clusterHits in pfamClusterHits.iteritems():
                
                line = ""
                
                line += str(pfamID)
               
                for clusterID in clusterHits.iterkeys():
                    
                    line += "\t" + str(clusterHits[clusterID])
                    
                line += "\t" + str(pfamClusterHitsMean[pfamID]) + "\t" + str(pfamClusterHitsVariance[pfamID]) + "\n"
                
                pfamClusterHitsStatsFile_fh.write(line)
                
            pfamClusterHitsStatsFile_fh.close()

        except Exception:
            
            print traceback.print_exc()


    def writeClusterPfamOrganismStatistics(self, clusterPfamOrganismHits, clusterPfamOrganismHitsFile):
        
        try:
            
            print "writeClusterPfamOrganismStatistics"
            
            clusterPfamOrganismHitsFile_fh = open(clusterPfamOrganismHitsFile, "w")
            
            header = "ClusterId"

            clusterPfamOrganismHitsKeys = clusterPfamOrganismHits.keys()
            
            pfamHitsPercentDict = clusterPfamOrganismHits[clusterPfamOrganismHitsKeys[0]][0]
            organismHitsPercentDict  = clusterPfamOrganismHits[clusterPfamOrganismHitsKeys[0]][1]

            for pfamID in pfamHitsPercentDict.iterkeys():
                
                header += "\t" + pfamID + "_Hits_%"
        
            for organismID in organismHitsPercentDict.iterkeys():
                
                header += "\t" + organismID + "_Hits_%"
                
                
            clusterPfamOrganismHitsFile_fh.write(header+"\n")

            line = ""
            
            for clusterID, clusterPfamOrgsHitsList in clusterPfamOrganismHits.iteritems():
                
                line = clusterID
                
                pfamHitsPercentDict = clusterPfamOrgsHitsList[0]
                organismHitsPercentDict  = clusterPfamOrgsHitsList[1]
                
                for pfamID in pfamHitsPercentDict.iterkeys():
                
                    line += "\t" + str(round(pfamHitsPercentDict[pfamID],5) )
                        
                for organismID in organismHitsPercentDict.iterkeys():
                
                    line += "\t" + str( round(organismHitsPercentDict[organismID],5) )
                    
                    
                clusterPfamOrganismHitsFile_fh.write(line+"\n")

                
            clusterPfamOrganismHitsFile_fh.close()

        except Exception:
            
            print traceback.print_exc()



###########################################################################----- Write to File ######################################################################################################



########################################################################### Pfam Cluster Organism Counts Final Output-------##########################################################################


#dict[pfam]= {cluster: {organisms_hit_dict{}, species_genes{}} }
#e.g. 
#dict[pfam] = { cluster1: [ {orgainsm[orgname1]:count,....., orgainsm[orgnam n]:count} , {orgainsm[orgnam 1]:[gene......gene],....., orgainsm[orgn n]:[gene,......,gene]}  ] }

    def createPfamClusterOrganismSelectedOrgsCountsDS(self, clusterIdsDict, pfamIdsDict):
        
        try:
            
            print "createPfamClusterOrganismSelectedOrgsCountsDS"
            
            pfamClusterOrganismCountsDict = {}
            
            for pfamID in pfamIdsDict.iterkeys():

                if not pfamClusterOrganismCountsDict.has_key(pfamID):
                    
                    pfamClusterOrganismCountsDict[pfamID] = {}
                
                clusterOrganismCounts = {}
                
                for clusterID in clusterIdsDict.iterkeys():
                    
                    if not clusterOrganismCounts.has_key(clusterID):
                            
                        clusterOrganismCounts[clusterID] = {}
                        
                    organismHitsCounts     = {}
                    selectedOrganismsGenes = {}
                    
                    for organismID in self.organismsSequenceIDsDict.iterkeys():
                    
                        if not organismHitsCounts.has_key(organismID):
                                
                            organismHitsCounts[organismID] = 0
                    
                    
                    for organismID in self.organismsToCompareList:
                    
                        if not selectedOrganismsGenes.has_key(organismID):
                                
                            selectedOrganismsGenes[organismID] = []
                            
                    
                    clusterOrganismCounts[clusterID]=[organismHitsCounts,selectedOrganismsGenes]
                    
                    
                pfamClusterOrganismCountsDict[pfamID] = clusterOrganismCounts
                
            return pfamClusterOrganismCountsDict
         
       
        
        except Exception:
            
            print traceback.print_exc()
            
        return ""
           
            
    def getPfamClusterOrganismSelectedOrgsStats(self, inflationValue, inflationClustersDict):
        
        try:
            
            print "getPfamClusterOrganismSelectedOrgsStats : " + str(inflationValue)
            
            pfamClusterOrganismCounts_inflation_Dict = {}

            pfamClusterOrganismCountsDict = {}
                        

            pfamIdsDict    = inflationClustersDict[inflationValue][0]
            clusterIdsDict = inflationClustersDict[inflationValue][1]
        
            pfamClusterOrganismCountsDict = self.createPfamClusterOrganismSelectedOrgsCountsDS(clusterIdsDict, pfamIdsDict)
                               
            for pfamID, pfamHits  in pfamIdsDict.iteritems():

                pfamClusters = pfamClusterOrganismCountsDict[pfamID]

                for clusterID, clusterHits in clusterIdsDict.iteritems():
                                            
                    clusterValuesList      = pfamClusters[clusterID]
                    organismHitsCounts     = clusterValuesList[0]
                    selectedOrganismsGenes = clusterValuesList[1]
                    
                    for protein in pfamHits:
                        
                        
                        if protein in clusterHits:
                                                    
                            for organism in organismHitsCounts.iterkeys():
                            
                                organismProteinsList = []
                                
                                organismProteinsList = self.organismsSequenceIDsDict[organism]
        
                                if protein in organismProteinsList:

                                    organismHitsCounts[organism] += 1
                                    
                                
                            for organism in self.organismsToCompareList:
                                
                                organismProteinsList = []
                                
                                organismProteinsList = self.organismsSequenceIDsDict[organism]
    
                                if protein in organismProteinsList:

                                    DataStrucures_Util.addUniqueListValueToDictionary(organism, protein, selectedOrganismsGenes)
                    
                                    
            pfamClusterOrganismCounts_inflation_Dict[inflationValue] = pfamClusterOrganismCountsDict

#            pfamClusterOrganismCountsFile = NGS_Util.createFilePath(self.clusterPfamOrganismStatsDir, "PfamClusterOrganismSelectedOrgsStats_" + inflationValue + ".txt")
            pfamClusterOrganismCountsFile = NGS_Util.createFilePath(self.interproResultDir, "PfamClusterOrganismSelectedOrgsStats_" + inflationValue + ".txt")
            self.writePfamClusterOrganismSelectedOrgsStats(pfamClusterOrganismCountsDict, pfamClusterOrganismCountsFile)

#            pfamClusterOrganismCountsFile = NGS_Util.createFilePath(self.clusterPfamOrganismStatsDir, "PfamClusterSelectedOrgsStats_" + inflationValue + ".txt")
            pfamClusterOrganismCountsFile = NGS_Util.createFilePath(self.interproResultDir, "PfamClusterSelectedOrgsStats_" + inflationValue + ".txt")
            self.writePfamClusterSelectedOrgsStats(pfamClusterOrganismCountsDict, pfamClusterOrganismCountsFile)
            
                
            return pfamClusterOrganismCounts_inflation_Dict

         
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def getNonInformativeInterprosPfamClusterOrganismSelectedOrgsStats(self, inflationValue, inflationClustersDict, interproResultDir):
        
        try:
            
            print "getNonInformativeInterprosPfamClusterOrganismSelectedOrgsStats : " + str(inflationValue)
            
            pfamClusterOrganismCounts_inflation_Dict = {}

            pfamClusterOrganismCountsDict = {}
                        

            pfamIdsDict    = inflationClustersDict[inflationValue][0]
            clusterIdsDict = inflationClustersDict[inflationValue][1]
        
            pfamClusterOrganismCountsDict = self.createPfamClusterOrganismSelectedOrgsCountsDS(clusterIdsDict, pfamIdsDict)
                               
            for pfamID, pfamHits  in pfamIdsDict.iteritems():

                pfamClusters = pfamClusterOrganismCountsDict[pfamID]

                for clusterID, clusterHits in clusterIdsDict.iteritems():
                                            
                    clusterValuesList      = pfamClusters[clusterID]
                    organismHitsCounts     = clusterValuesList[0]
                    selectedOrganismsGenes = clusterValuesList[1]
                    
                    for protein in pfamHits:
                        
                        
                        if protein in clusterHits:
                                                    
                            for organism in organismHitsCounts.iterkeys():
                            
                                organismProteinsList = []
                                
                                organismProteinsList = self.organismsSequenceIDsDict[organism]
        
                                if protein in organismProteinsList:

                                    organismHitsCounts[organism] += 1
                                    
                                
                            for organism in self.organismsToCompareList:
                                
                                organismProteinsList = []
                                
                                organismProteinsList = self.organismsSequenceIDsDict[organism]
    
                                if protein in organismProteinsList:

                                    DataStrucures_Util.addUniqueListValueToDictionary(organism, protein, selectedOrganismsGenes)
                    
                                    
            pfamClusterOrganismCounts_inflation_Dict[inflationValue] = pfamClusterOrganismCountsDict

#            pfamClusterOrganismCountsFile = NGS_Util.createFilePath(self.clusterPfamOrganismStatsDir, "PfamClusterOrganismSelectedOrgsStats_" + inflationValue + ".txt")
            pfamClusterOrganismCountsFile = NGS_Util.createFilePath(interproResultDir, "PfamClusterOrganismSelectedOrgsStats_" + inflationValue + ".txt")
            self.writePfamClusterOrganismSelectedOrgsStats(pfamClusterOrganismCountsDict, pfamClusterOrganismCountsFile)

#            pfamClusterOrganismCountsFile = NGS_Util.createFilePath(self.clusterPfamOrganismStatsDir, "PfamClusterSelectedOrgsStats_" + inflationValue + ".txt")
            pfamClusterOrganismCountsFile = NGS_Util.createFilePath(interproResultDir, "PfamClusterSelectedOrgsStats_" + inflationValue + ".txt")
            self.writePfamClusterSelectedOrgsStats(pfamClusterOrganismCountsDict, pfamClusterOrganismCountsFile)
            
                
            return pfamClusterOrganismCounts_inflation_Dict

         
        except Exception:
            
            print traceback.print_exc()
            
        return ""



    def writePfamClusterOrganismSelectedOrgsStats(self, pfamClusterOrganismCountsDict, pfamClusterOrganismCountsFile):
        
        try:
            
            print "writePfamClusterOrganismSelectedOrgsStats"
            
            pfamClusterOrganismSelectedOrgsStats_fh = open(pfamClusterOrganismCountsFile, "w")

            organismsSequenceIDsList =  self.organismsSequenceIDsDict.keys()
            
            header = "PfamID\tClusterId"

            for organismID in organismsSequenceIDsList:
                
                header += "\t" + organismID
                
            for organismID in self.organismsToCompareList:
                
                header += "\t" + organismID
        
                
            pfamClusterOrganismSelectedOrgsStats_fh.write(header+"\n")

            line = ""
            
            for pfamID, pfamClusterOrganismCountsValues in pfamClusterOrganismCountsDict.iteritems():
                              
                for clusterID, clusterOrganismCountsValues in pfamClusterOrganismCountsValues.iteritems():
                    
                    line = pfamID + "\t" + clusterID
                    
                    organismHitsCounts     = clusterOrganismCountsValues[0]
                    selectedOrganismsGenes = clusterOrganismCountsValues[1]
                    
                    for organismName in organismsSequenceIDsList:
                        
                        line += "\t" + str(organismHitsCounts[organismName])
                                                    
                    for organismName in self.organismsToCompareList:
                        
                        line += "\t" + ",".join(selectedOrganismsGenes[organismName])
                    
                    pfamClusterOrganismSelectedOrgsStats_fh.write(line+"\n")

                
            pfamClusterOrganismSelectedOrgsStats_fh.close()

        except Exception:
            
            print traceback.print_exc()
        


    def writePfamClusterSelectedOrgsStats(self, pfamClusterOrganismCountsDict, pfamClusterOrganismCountsFile):
        
        try:
            
            print "writePfamClusterSelectedOrgsStats"
            
            pfamClusterOrganismSelectedOrgsStats_fh = open(pfamClusterOrganismCountsFile, "w")

            
            header = "PfamID\tClusterId"

            for organismID in self.organismsToCompareList:
                
                header += "\t" + organismID
        
                
            pfamClusterOrganismSelectedOrgsStats_fh.write(header+"\n")

            line = ""
            
            for pfamID, pfamClusterOrganismCountsValues in pfamClusterOrganismCountsDict.iteritems():
                              
                for clusterID, clusterOrganismCountsValues in pfamClusterOrganismCountsValues.iteritems():
                    
                    line = pfamID + "\t" + clusterID
                    
                    organismHitsCounts     = clusterOrganismCountsValues[0]
                    selectedOrganismsGenes = clusterOrganismCountsValues[1]
                    
                    for organismName in self.organismsToCompareList:
                        
                        line += "\t" + ",".join(selectedOrganismsGenes[organismName])
                    
                    pfamClusterOrganismSelectedOrgsStats_fh.write(line+"\n")

                
            pfamClusterOrganismSelectedOrgsStats_fh.close()

        except Exception:
            
            print traceback.print_exc()
        
###########################################################################----- Pfam Cluster Organism Counts Final Output############################################################################




########################################################################### Interpro Pfam Organism Counts Final Output-------##########################################################################


#dict[pfam]= {cluster: {organisms_hit_dict{}, species_genes{}} }
#e.g. 
#dict[pfam] = { cluster1: [ {orgainsm[orgname1]:count,....., orgainsm[orgnam n]:count} , {orgainsm[orgnam 1]:[gene......gene],....., orgainsm[orgn n]:[gene,......,gene]}  ] }

    def createPfamOrganismSelectedOrgsCountsDS(self, pfamIdsDict):
        
        try:
            
            print "createPfamOrganismSelectedOrgsCountsDS"
            
            pfamOrganismCountsDict = {}
            
            for pfamID in pfamIdsDict.iterkeys():

                if not pfamOrganismCountsDict.has_key(pfamID):
                    
                    pfamOrganismCountsDict[pfamID] = []
                                        
                    organismHitsCounts     = {}
                    selectedOrganismsGenes = {}
                    
                    for organismID in self.organismsSequenceIDsDict.iterkeys():
                    
                        if not organismHitsCounts.has_key(organismID):
                                
                            organismHitsCounts[organismID] = 0
                    
                    
                    for organismID in self.organismsToCompareList:
                    
                        if not selectedOrganismsGenes.has_key(organismID):
                                
                            selectedOrganismsGenes[organismID] = []
                                
                                           
                    pfamOrganismCountsDict[pfamID] = [organismHitsCounts,selectedOrganismsGenes]
                
            return pfamOrganismCountsDict
         
       
        
        except Exception:
            
            print traceback.print_exc()
            
        return ""
           
############################## to work below ## fix it             
    def getInterproPfamOrganismSelectedOrgsStats(self, pfamID, bestPfamInflationClusterValues, bestPfamClusterOrganismCountsDict):
        
        try:
            
            print "getInterproPfamOrganismSelectedOrgsStats"
            

            pfamOrganismStatsDict = {}
            interproCorrespondingPfamDict = {}

            pfamIdsDict                           =  bestPfamInflationClusterValues[0]
            clusterIdsDict                        =  bestPfamInflationClusterValues[1]            
            interproCorrespondingPfamDict[pfamID] =  pfamIdsDict[pfamID]           
            pfamOrganismStatsDict                 =  self.createPfamOrganismSelectedOrgsCountsDS(interproCorrespondingPfamDict)
            
            pfamClusterOrganismCountsValues       =  bestPfamClusterOrganismCountsDict[pfamID]
            
            
#            for pfamID, pfamClusterOrganismCountsValues in bestPfamClusterOrganismCountsDict.iteritems():
                
            pfamOrganismStatsValues = pfamOrganismStatsDict[pfamID]
                          
            for clusterID, clusterOrganismCountsValues in pfamClusterOrganismCountsValues.iteritems():
                
                line = pfamID + "\t" + clusterID
                
                bestPfamOrganismHitsCounts     = clusterOrganismCountsValues[0]
                bestPfamSelectedOrganismsGenes = clusterOrganismCountsValues[1]
                
                for organismName in self.organismsSequenceIDsDict.iterkeys():
                    
                    pfamOrganismStatsValues[0][organismName] += bestPfamOrganismHitsCounts[organismName]
                                                
                for organismName in self.organismsToCompareList:
                    
                    pfamOrganismStatsValues[1][organismName] = DataStrucures_Util.mergeListsToUniqueValuesList(pfamOrganismStatsValues[1][organismName],bestPfamSelectedOrganismsGenes[organismName])
                
                            
            return pfamOrganismStatsDict

         
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def writeInterproPfamOrganismSelectedOrgsStats(self):
        
        try:
            
            print "writeInterproPfamOrganismSelectedOrgsStats"
            
            interproPfamOrganismSelectedOrgsStatsFile = NGS_Util.createFilePath(self.resultDir, "InterproPfamOrganismStats.txt")

            interproPfamOrganismSelectedOrgsStatsFile_fh = open(interproPfamOrganismSelectedOrgsStatsFile, "w")
           
            header = "InterproID\tInterproDescription\tPfamId"

            for organismID in self.organismsSequenceIDsDict.keys():
                
                header += "\t" + organismID
                
            for organismID in self.organismsToCompareList:
                
                header += "\t" + organismID
                        
            interproPfamOrganismSelectedOrgsStatsFile_fh.write(header+"\n")


            line = ""
            
            
######            self.interproPfamOrganismSelectedOrgsStatsDict[interproID] = self.getInterproPfamOrganismSelectedOrgsStats(fullPfamID, bestPfamInflationClusterValues, bestPfamClusterOrganismCountsDict)
            
            for interproID, pfamOrganismStatsDict in self.interproPfamOrganismSelectedOrgsStatsDict.iteritems():
                
                if "_" in interproID:
                    interproID = interproID[0:interproID.find("_",0)]
                    
                
                for pfamID, pfamOrganismStatsValues in pfamOrganismStatsDict.iteritems():
                    
                    line = interproID + "\t" + self.interproDescDict[interproID]  +  "\t"+ pfamID
                    
                    organismHitsCounts     = pfamOrganismStatsValues[0]
                    selectedOrganismsGenes = pfamOrganismStatsValues[1]
                    
                    for organismName in self.organismsSequenceIDsDict.keys():
                        
                        line += "\t" + str(organismHitsCounts[organismName])
                                                    
                    for organismName in self.organismsToCompareList:
                        
                        line += "\t" + ",".join(selectedOrganismsGenes[organismName])
                    
                    interproPfamOrganismSelectedOrgsStatsFile_fh.write(line+"\n")

                
            interproPfamOrganismSelectedOrgsStatsFile_fh.close()

        except Exception:
            
            print traceback.print_exc()


###########################################################################----- Interpro Pfam Organism Counts Final Output############################################################################


    def binAverageInflationValue(self, averageInflationValue):
        
        try:
            
            print "binAverageInflationValue + "
            

            inflationValue = 1.2
            
            binnedInflationValue = 1.2
            
            averageInflationValue = float(averageInflationValue)
            
            minDiff = 0.4
            
            for index in range(1,10):
                
                diff = abs(averageInflationValue - inflationValue)
                
                if  diff == 0:
                    
                    binnedInflationValue = inflationValue
                    break
                
                else:
                    if diff < minDiff:
                        
                        minDiff = diff
                        binnedInflationValue = inflationValue
                        
                
                inflationValue += 0.4
                

            return binnedInflationValue

            
        except Exception:
            
            print traceback.print_exc()
        
            return ""




    def getNonInformativePfamsAverageInflationValue(self):
        
        try:
            
            print "getNonInformativePfamsAverageInflationValue"
            
            inflationValue = 0.0

            for interproID, interproValuesList  in self.interproInflationValuesSelectedDict.iteritems():
                               
                inflationValue += (float(self.interproInflationValuesSelectedDict[interproID][1])*1.0)/10
                
            
            if inflationValue > 0.0:                         
                
                averageInflationValue = inflationValue/len(self.interproInflationValuesSelectedDict)
            else:
                averageInflationValue = 0.0
                
            return averageInflationValue

            
        except Exception:
            
            print traceback.print_exc()
        
            return ""


    def calculateNonInformativePfamsInflationValue(self):
    
        try:

            print "calculateNonInformativePfamsInflationValue"
            
            averageInflationValue = self.getNonInformativePfamsAverageInflationValue()
            
            binnedInflationValue = str(self.binAverageInflationValue(averageInflationValue)).replace(".","")
            
            if averageInflationValue > 0.0:
            
                for interproID, interproValuesList  in self.nonInformativeInterproDict.iteritems():
    
                    pfamStatsFile = NGS_Util.createFilePath(interproValuesList[1], "average_pfam_inflation_value.txt")
                    pfamStatsFile_fh = open(pfamStatsFile, "w")
                    
                    pfamStatsFile_fh.write( "InterproID\tAverageBinnedInflationValue\n")
                    pfamStatsFile_fh.write( interproID + "\t" + binnedInflationValue + "\n")
                        
                    pfamStatsFile_fh.close()
                    
                    
                    inflationClustersDict = interproValuesList[3]
                                   
                    pfamClusterOrganismCounts_inflation_Dict = self.getNonInformativeInterprosPfamClusterOrganismSelectedOrgsStats(binnedInflationValue, inflationClustersDict, interproValuesList[1])
                    
                    #pfamClusterOrganismCounts_inflation_Dict = interproValuesList[3]
                    
                    inflationValue = binnedInflationValue #str(binnedInflationValue).replace(".","")
                    
                    bestPfamInflationClusterValues    = inflationClustersDict[inflationValue]
                    bestPfamClusterOrganismCountsDict = pfamClusterOrganismCounts_inflation_Dict[inflationValue]
                    
                    self.interproPfamOrganismSelectedOrgsStatsDict[interproID] = self.getInterproPfamOrganismSelectedOrgsStats(interproValuesList[0], bestPfamInflationClusterValues, bestPfamClusterOrganismCountsDict)
            
            else:
                
                for interproID, interproValuesList  in self.nonInformativeInterproDict.iteritems():
                
                    fullPfamID            = interproValuesList[0]
                    inflationClustersDict = interproValuesList[3]
                    mostVariedPfamList    = interproValuesList[4]
                
                    mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff = self.doInflationValueCalculations(inflationClustersDict, mostVariedPfamList) 
                    self.interproInflationValuesSelectedDict[interproID] = [mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff]
    
                    pfamClusterOrganismCounts_inflation_Dict = self.getNonInformativeInterprosPfamClusterOrganismSelectedOrgsStats(binnedInflationValue, inflationClustersDict, interproValuesList[1])
                
                    bestPfamInflationClusterValues    = inflationClustersDict[bestInflationValue]
                    bestPfamClusterOrganismCountsDict = pfamClusterOrganismCounts_inflation_Dict[bestInflationValue]
                    
                    self.interproPfamOrganismSelectedOrgsStatsDict[interproID] = self.getInterproPfamOrganismSelectedOrgsStats(fullPfamID, bestPfamInflationClusterValues, bestPfamClusterOrganismCountsDict)
                    
                
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def getInterproPfamsMAD(self, pfamID, blastResultFile):
        
        try:
            
            print "doSpecificitySensitivityCalculations"
            
            pfamBlastStatistics = PfamBlastStatistics.PfamBlastStatistics()
            pfamBlastStatistics.initialize(self.resultDir, self.interproResultDir, self.blastResultDir, self.pfamBlastStaticticsDir)
            
            pfamMed, pfamMAD = pfamBlastStatistics.getPfamMedianAbsoluteDeviation( pfamID, blastResultFile)
            
            
            return pfamMed, pfamMAD

        except Exception:
            
            print traceback.print_exc()
            
            
    def insertInterproPfamMADDict(self, interproID, pfamID,pfamMed, pfamMAD):
        
        try:
            
            print "doSpecificitySensitivityCalculations"
            
            valuesDict = {}
            
            
            if not self.interproPfamMADDict.has_key(interproID):
                
                self.interproPfamMADDict[interproID] = valuesDict

            
            valuesDict = self.interproPfamMADDict[interproID]
            
            valuesDict[pfamID] = [pfamMed, pfamMAD]
                


        except Exception:
            
            print traceback.print_exc()
            

    def writeInterproPfamsMAD(self):
        
        try:
            
            print "doSpecificitySensitivityCalculations"
            
            interproPfamMADFile = NGS_Util.createFilePath(self.resultDir, "InterproPfamMAD.txt")

            header = "InterproID\tPfamId\tMedian\tMeanAbsoluteDeviation"
            
            interproPfamMADFile_fh = open(interproPfamMADFile, "w")
            
            valuesDict = {}
            
            line = ""
            
            interproPfamMADFile_fh.write(header + "\n")
            
            for interproID, valuesDict in self.interproPfamMADDict.iteritems():
                               
                for pfamID, pfamValues in valuesDict.iteritems():
                    
                    line = interproID + "\t" + pfamID
                    
                    for val in pfamValues:
                        
                        line += "\t" + str(val)
                        
                    interproPfamMADFile_fh.write(line + "\n")
                
            interproPfamMADFile_fh.close()

        except Exception:
            
            print traceback.print_exc()
                        


################################################################################ compare pfam files-----###################################################################################

# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>

#FGSG.03090     46    206     46    206 PF01554.13  MatE              Family     1   162   162    101.7     3e-29   1 CL0222   

    def getPfamScanHitsList(self, pfamScan_output, pfamScanHitsList):
    
        try:

#            print "getPfamScanHitsList: " + pfamScan_output
            
            pfamScan_output_fh = open(pfamScan_output)
            
            line = ""
            
            hitsList=[]
            
            for line in pfamScan_output_fh:

                line = line.strip()

                if not line.startswith("#"):
                    
                    if len(line) > 0:
                        
                        hitsList = line.split(" ")
                                       
                        hitsList = [val for val in hitsList if len(val) > 0]
                        
                        pfamScanHitsList.append(hitsList)
                                                            
            pfamScan_output_fh.close()
            
            return pfamScanHitsList
     
        except Exception:
            
            print traceback.print_exc()
            
        return pfamScanHitsList


    def getMockPfamScanHits(self,pfamID):
    
        try:

            print "getMockPfamScanHits"  + pfamID
            
            mockPfamScanHitsList = []
                    
            orgListFile_fh = open(self.orgListFile)

            for orgLine in orgListFile_fh:
                
                organismName = orgLine.strip()
                
                org_pfamScan = NGS_Util.createFilePath(self.orgPfamResultDir, organismName + "." + pfamID + ".pfamScan")                
                
                mockPfamScanHitsList = self.getPfamScanHitsList(org_pfamScan, mockPfamScanHitsList)

            orgListFile_fh.close()
            
            
            print mockPfamScanHitsList
            
            
            return mockPfamScanHitsList
     
        except Exception:
            
            print traceback.print_exc()
            
        return []


    def getWholePfamScanHits(self,pfamID):
    
        try:

            print "getWholePfamScanHits:" + pfamID
            
            pfamScanHitsList = []

            pfamScan = NGS_Util.createFilePath(self.interproResultDir, pfamID + ".pfamScan")
                                    
            pfamScanHitsList = self.getPfamScanHitsList(pfamScan, pfamScanHitsList)
            
            print pfamScanHitsList
            
            return pfamScanHitsList
     
        except Exception:
            
            print traceback.print_exc()
            
        return []


    def findMatchingPfaHits(self,pfamID):
    
        try:

            print "findMatchingPfaHits"

            pfamScanHitsList     = []            
            mockPfamScanHitsList = []
            matchingPfamHits     = []

            
            pfamScanHitsList     = self.getWholePfamScanHits(pfamID)
            mockPfamScanHitsList = self.getMockPfamScanHits(pfamID)
            

            for pfamScanHit in pfamScanHitsList:
            
                for mockPfamScanHit in mockPfamScanHitsList:
                    
                    print "comparing pfam: " + pfamScanHit[5] + " ..." + mockPfamScanHit[5]
                
                    if pfamScanHit[5] == mockPfamScanHit[5]:
                        
                        matchingPfamHits.append(pfamScanHit)
                
            
            return matchingPfamHits
            
        except Exception:
            
            print traceback.print_exc()
            
        return []

                
    def writeMatchingPfamscanHits(self, pfamID,matchingPfamHits, matchingPfamHitsFile ):
    
        
        try:
            
            print "writeMatchingPfamscanHits"
            
            pfamScan_output = NGS_Util.createFilePath(self.interproResultDir, pfamID + ".pfamScan")
            
            pfamScan_output_fh = open(pfamScan_output)
            
            matchingPfamHitsFile_fh = open(matchingPfamHitsFile, "w")
            
            
            line = ""
            
           
            for line in pfamScan_output_fh:

                if line.startswith("#"):
                    
                    matchingPfamHitsFile_fh.write(line)
            
            matchingPfamHitsFile_fh.write("\n")

                    
            for hit in matchingPfamHits:
                
                for val in hit:
                    
                    line += str(val) + "\t"
                    
                
                matchingPfamHitsFile_fh.write(line.strip() + "\n")
            
                
            matchingPfamHitsFile_fh.close()
            
            pfamScan_output_fh.close()
                    
        except Exception:
            
            print traceback.print_exc()



    def comparePfamsScanResults(self,pfamID):
    
        try:
            
            print "removeInterproScanWithUnmatchedPfamsScanResults"

            matchingPfamHits = self.findMatchingPfaHits(pfamID)
            
            if len(matchingPfamHits) > 0:
            
                return True
            
            else:
                
                print "deleting dir:" +  self.interproResultDir
                shutil.rmtree(self.interproResultDir)
                #os.removedirs(self.interproResultDir)
            
            return False
     
        except Exception:
            
            print traceback.print_exc()


    #
    #def comparePfamsScanResults(self,pfamID):
    #
    #    try:
    #        
    #        print "removeInterproScanWithUnmatchedPfamsScanResults"
    #
    #        matchingPfamHitsFile = NGS_Util.createFilePath(self.interproResultDir, pfamID + "_corrected.pfamScan")
    #        
    #        matchingPfamHits = self.findMatchingPfaHits(pfamID)
    #        
    #        if len(matchingPfamHits) > 0:
    #        
    #            self.writeMatchingPfamscanHits(pfamID,matchingPfamHits, matchingPfamHitsFile )
    #            return matchingPfamHitsFile
    #        
    #        else:
    #            
    #            print "deleting dir:" +  self.interproResultDir
    #            #os.removedirs(self.interproResultDir)
    #        
    #        return ""
    # 
    #    except Exception:
    #        
    #        print traceback.print_exc()
################################################################################----- compare pfam files###################################################################################


    def runForInterPro(self, interproID):
    
        try:

            print "runAll"
            
            
            if self.interproPfamsDict.has_key(interproID):
                
                pfamIDsList =  self.interproPfamsDict[interproID]
                                
                for pfamID in pfamIDsList:
                    
                    if len(pfamIDsList) > 1 :
                        self.makeInterproResultsDirectories(interproID + "_" + pfamID, self.resultDir)
                        interproRunID = interproID + "_" + pfamID
                    else:
                        self.makeInterproResultsDirectories(interproID, self.resultDir)
                        interproRunID = interproID 
                        
                    print "Analysing organisms for interpro id: " + interproID + " and PfamID " +  pfamID
                    
                    
                    fullPfamID = self.grepPfamID_from_PfamA_hmm(pfamID,self.PfamA_hmm_file)
                    
                    self.createMockPfam(pfamID, fullPfamID, self.PfamA_hmm_file)            
                    
                    #self.doOrgPfamScan(pfamID)
                    self.doThreadedOrgPfamScan(pfamID)
                    
                    pfamSearchResultSequenceFile = NGS_Util.createFilePath(self.interproResultDir, pfamID + ".faa")
                    pfamSearchResultSequenceFile_fh = open(pfamSearchResultSequenceFile,"w")
                    pfamSearchResultSequenceFile_fh.close()
                    
                    
                    self.makeBlast_PfamScanSequence_File(pfamID, pfamSearchResultSequenceFile)
        
                    blastResultFile = self.doPfamScanSequenceBlast(pfamID, pfamSearchResultSequenceFile)
                    
                    threshold = self.getThreshold(pfamID, blastResultFile)
                    
                    eValueCutOff = self.getEValueCutOff(threshold)
                    
                    mclInputFile = self.makeMclInputFile(pfamID, blastResultFile, eValueCutOff)
                    
                    self.doMCLClustering(pfamID, mclInputFile)
        
                    pfamScan = self.doPfamScan(pfamID, pfamSearchResultSequenceFile)
                    
                    
                    matchingPfamHitsFile = self.comparePfamsScanResults(pfamID)
                    
                    
                    if matchingPfamHitsFile:                    
                    
                        inflationClustersDict, mostVariedPfamList, pfamVarianceDict = self.doSpecificitySensitivityCalculations(pfamScan, pfamID, fullPfamID)
                                            
                        pfamIdsDict = inflationClustersDict.values()[0][0]
                                   
                        if len(pfamIdsDict) > 1:
                            
                            mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff = self.doInflationValueCalculations(inflationClustersDict, mostVariedPfamList) 
                            self.interproInflationValuesSelectedDict[interproRunID] = [mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff]
            
                            pfamClusterOrganismCounts_inflation_Dict = self.getPfamClusterOrganismSelectedOrgsStats(bestInflationValue, inflationClustersDict)
                        
                            bestPfamInflationClusterValues    = inflationClustersDict[bestInflationValue]
                            bestPfamClusterOrganismCountsDict = pfamClusterOrganismCounts_inflation_Dict[bestInflationValue]
                            
                            self.interproPfamOrganismSelectedOrgsStatsDict[interproRunID] = self.getInterproPfamOrganismSelectedOrgsStats(fullPfamID, bestPfamInflationClusterValues, bestPfamClusterOrganismCountsDict)
                            
                            pfamMed, pfamMAD = self.getInterproPfamsMAD(pfamID, blastResultFile)
                            self.insertInterproPfamMADDict(interproID, pfamID,pfamMed, pfamMAD)
                                                    
                                    
                        elif len(pfamIdsDict) == 1:
                            
                            pfamMed, pfamMAD = self.getInterproPfamsMAD(pfamID, blastResultFile)
                            self.insertInterproPfamMADDict(interproID, pfamID,pfamMed, pfamMAD)
    
                            self.nonInformativeInterproDict[interproRunID] = [fullPfamID, self.interproResultDir, self.SensitivitySpecificityDir, inflationClustersDict, mostVariedPfamList]
                    

        except Exception:
            
            print traceback.print_exc()
            
        return ""




    def runForAllInterPros(self):
    
        try:

            print "runForAllInterPros"
                        
            for interproID in self.interproPfamsDict.iterkeys():
                
                self.runForInterPro(interproID)

        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def redoInterproCalculations(self, interproID):
    
        try:

            print "redoInterproCalculations"
            
            
            if self.interproPfamsDict.has_key(interproID):
                
                pfamIDsList =  self.interproPfamsDict[interproID]
                                
                for pfamID in pfamIDsList:
                    
                    if len(pfamIDsList) > 1 :
                        self.makeInterproResultsDirectoriesPath(interproID + "_" + pfamID, self.resultDir)
                        interproRunID = interproID + "_" + pfamID
                    else:
                        self.makeInterproResultsDirectoriesPath(interproID, self.resultDir)
                        interproRunID = interproID 
                        
                    if os.path.exists(self.interproResultDir):
                        
                        print "Analysing organisms for interpro id: " + interproID + " and PfamID " +  pfamID
                       
                        
                        fullPfamID = self.grepPfamID_from_PfamA_hmm(pfamID,self.PfamA_hmm_file)
                        
                        pfamSearchResultSequenceFile = NGS_Util.createFilePath(self.interproResultDir, pfamID + ".faa")
                        
                        blastResultFile = NGS_Util.createFilePath(self.blastResultDir, pfamID + ".blast")
                        
                        pfamScan = NGS_Util.createFilePath(self.interproResultDir, pfamID + ".pfamScan")
                        
                        matchingPfamHitsFile = self.comparePfamsScanResults(pfamID)
                        
                        
                        if matchingPfamHitsFile:
                        
                            inflationClustersDict, mostVariedPfamList, pfamVarianceDict = self.doSpecificitySensitivityCalculations(pfamScan, pfamID, fullPfamID)
                                                
                            pfamIdsDict = inflationClustersDict.values()[0][0]
                                       
                            if len(pfamIdsDict) > 1:
                                
                                mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff = self.doInflationValueCalculations(inflationClustersDict, mostVariedPfamList) 
                                self.interproInflationValuesSelectedDict[interproRunID] = [mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff]
                            
                                pfamClusterOrganismCounts_inflation_Dict = self.getPfamClusterOrganismSelectedOrgsStats(bestInflationValue, inflationClustersDict)
                            
                                bestPfamInflationClusterValues    = inflationClustersDict[bestInflationValue]
                                bestPfamClusterOrganismCountsDict = pfamClusterOrganismCounts_inflation_Dict[bestInflationValue]
                                
                                self.interproPfamOrganismSelectedOrgsStatsDict[interproRunID] = self.getInterproPfamOrganismSelectedOrgsStats(fullPfamID, bestPfamInflationClusterValues, bestPfamClusterOrganismCountsDict)
                                
                                pfamMed, pfamMAD = self.getInterproPfamsMAD(pfamID, blastResultFile)
                                self.insertInterproPfamMADDict(interproID, pfamID,pfamMed, pfamMAD)
                                                        
                                        
                            elif len(pfamIdsDict) == 1:
                                
                                pfamMed, pfamMAD = self.getInterproPfamsMAD(pfamID, blastResultFile)
                                self.insertInterproPfamMADDict(interproID, pfamID,pfamMed, pfamMAD)
                            
                                self.nonInformativeInterproDict[interproRunID] = [fullPfamID, self.interproResultDir, self.SensitivitySpecificityDir, inflationClustersDict, mostVariedPfamList]
                    
                    
        except Exception:
            
            print traceback.print_exc()
            
        return ""



    #def redoInterproCalculations(self, interproID):
    #
    #    try:
    #
    #        print "redoInterproCalculations"
    #        
    #        
    #        if self.interproPfamsDict.has_key(interproID):
    #            
    #            pfamIDsList =  self.interproPfamsDict[interproID]
    #                            
    #            for pfamID in pfamIDsList:
    #                
    #                if len(pfamIDsList) > 1 :
    #                    self.makeInterproResultsDirectories(interproID + "_" + pfamID, self.resultDir)
    #                    interproRunID = interproID + "_" + pfamID
    #                else:
    #                    self.makeInterproResultsDirectories(interproID, self.resultDir)
    #                    interproRunID = interproID 
    #                    
    #                print "Analysing organisms for interpro id: " + interproID + " and PfamID " +  pfamID
    #                
    #                
    #                fullPfamID = self.grepPfamID_from_PfamA_hmm(pfamID,self.PfamA_hmm_file)
    #                
    #                print "---------------fulll pfam id" + fullPfamID
    #                
    #                pfamSearchResultSequenceFile = NGS_Util.createFilePath(self.interproResultDir, pfamID + ".faa")
    #                
    #                blastResultFile = NGS_Util.createFilePath(self.blastResultDir, pfamID + ".blast")
    #                
    #                pfamScan = NGS_Util.createFilePath(self.interproResultDir, pfamID + ".pfamScan")
    #                
    #                inflationClustersDict, mostVariedPfamList, pfamVarianceDict = self.doSpecificitySensitivityCalculations(pfamScan, pfamID, fullPfamID)
    #                                    
    #                pfamIdsDict = inflationClustersDict.values()[0][0]
    #                           
    #                if len(pfamIdsDict) > 1:
    #                    
    #                    mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff = self.doInflationValueCalculations(inflationClustersDict, mostVariedPfamList) 
    #                    self.interproInflationValuesSelectedDict[interproRunID] = [mostVariedPfamID, bestInflationValue, min_specficity_sensitivity_diff]
    #                
    #                    pfamClusterOrganismCounts_inflation_Dict = self.getPfamClusterOrganismSelectedOrgsStats(bestInflationValue, inflationClustersDict)
    #                
    #                    bestPfamInflationClusterValues    = inflationClustersDict[bestInflationValue]
    #                    bestPfamClusterOrganismCountsDict = pfamClusterOrganismCounts_inflation_Dict[bestInflationValue]
    #                    
    #                    self.interproPfamOrganismSelectedOrgsStatsDict[interproRunID] = self.getInterproPfamOrganismSelectedOrgsStats(fullPfamID, bestPfamInflationClusterValues, bestPfamClusterOrganismCountsDict)
    #                    
    #                    pfamMed, pfamMAD = self.getInterproPfamsMAD(pfamID, blastResultFile)
    #                    self.insertInterproPfamMADDict(interproID, pfamID,pfamMed, pfamMAD)
    #                                            
    #                            
    #                elif len(pfamIdsDict) == 1:
    #                    
    #                    pfamMed, pfamMAD = self.getInterproPfamsMAD(pfamID, blastResultFile)
    #                    self.insertInterproPfamMADDict(interproID, pfamID,pfamMed, pfamMAD)
    #                
    #                    self.nonInformativeInterproDict[interproRunID] = [fullPfamID, self.interproResultDir, self.SensitivitySpecificityDir, inflationClustersDict, mostVariedPfamList]
    #                
    #
    #    except Exception:
    #        
    #        print traceback.print_exc()
    #        
    #    return ""


    def redoAllInterProsCalculations(self):
    
        try:

            print "redoAllInterProsCalculations"
                        
            for interproID in self.interproPfamsDict.iterkeys():
                
                self.redoInterproCalculations(interproID)

        except Exception:
            
            print traceback.print_exc()
            
        return ""



    def resumeRunForAllInterPros(self):
   
       try:

           print "resumeRunForAllInterPros"
                       
           for interproID in self.interproPfamsDict.iterkeys():
               
               iPath = NGS_Util.createDirectoryPath(self.resultDir,interproID)
               
               if not os.path.exists(iPath):
                   
                   self.runForInterPro(interproID)
           
           self.redoAllInterProsCalculations()

       except Exception:
            
            print traceback.print_exc()
            
       return ""



#    def resumeRunForAllInterPros(self):
    
#        try:

#            print "resumeRunForAllInterPros"
                        
#            for interproID in self.interproPfamsDict.iterkeys():
                
#                pfamIDsList =  self.interproPfamsDict[interproID]
                
#                for pfamID in pfamIDsList:
                    
#                    if len(pfamIDsList) > 1 :
#                        interproRunID = interproID + "_" + pfamID
#                    else:
#                        interproRunID = interproID
#                        
#                    iPath = NGS_Util.createDirectoryPath(self.resultDir,interproRunID)
                
#                    if not os.path.exists(iPath):
                    
#                        self.runForInterPro(interproID)
            
#            self.redoAllInterProsCalculations()

#        except Exception:
            
#            print traceback.print_exc()
            
#        return ""


