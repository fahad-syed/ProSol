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

class PfamBlastStatistics:

    resultDir              = ""
    interproResultDir      = ""
    blastResultDir         = ""
    pfamBlastStaticticsDir = ""
    
    interproPfamDivergenceScores = ""
    
            
    #Blast tabular output fields
    #0:  query name
    #1:  subject name
    #2:  percent identities
    #3:  aligned length
    #4:  number of mismatched positions
    #5:  number of gap positions
    #6:  query sequence start
    #7:  query sequence end
    #8:  subject sequence start
    #9:  subject sequence end
    #10: e-value
    #11: bit score

    def initialize(self, resultDir, interproResultDir, blastResultDir, pfamBlastStaticticsDir):
    
        try:
            
            self.resultDir              = resultDir
            self.interproResultDir      = interproResultDir
            self.blastResultDir         = blastResultDir
            self.pfamBlastStaticticsDir = pfamBlastStaticticsDir
            
        except Exception:
            
            print traceback.print_exc()


    def createBlastBitScoreMatrix(self, blastResultFile):
    
        try:

            print "createBlastBitScoreMatrix"
            
            bitScoreMatrix = {}
            
            if os.path.exists(blastResultFile):
                
                blastResultFile_fh = open(blastResultFile)
                
                line = ""
                blastOutput = []

       
                for line in blastResultFile_fh:
                    
                    blastOutput = line.strip().split("\t")                    
                    sequenceID  = blastOutput[0]
                    
                    if len(sequenceID) > 0:
                        
                        if not bitScoreMatrix.has_key(sequenceID):
                        
                            bitScoreMatrix[sequenceID]={}


                matchingSequenceIDsList = sorted(bitScoreMatrix.keys())

                
                for sequenceId in bitScoreMatrix.iterkeys():
            
                    for matchingSequenceID in matchingSequenceIDsList:

                        if not bitScoreMatrix[sequenceId].has_key(matchingSequenceID):
                        
                            bitScoreMatrix[sequenceId][matchingSequenceID] = 0.0

    
                blastResultFile_fh.close()
                
            return bitScoreMatrix
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def getBlastBitScoreMatrix(self, pfamID, blastResultFile):
    
        try:

            print "getBlastBitScoreMatrix"
            
            bitScoreMatrix = self.createBlastBitScoreMatrix(blastResultFile)
            
            if os.path.exists(blastResultFile):
                
                blastResultFile_fh = open(blastResultFile)
   
                line = ""            
                blastOutput = []
    
                for line in blastResultFile_fh:
                    
                    blastOutput = line.strip().split("\t")
                    
                    if len(blastOutput)>1:
    
                        sequenceID         = blastOutput[0]
                        matchingSequenceID = blastOutput[1]
                        bitScore           = float(blastOutput[11])
                    
                        if bitScoreMatrix[sequenceID][matchingSequenceID] < bitScore:
                            bitScoreMatrix[sequenceID][matchingSequenceID] = bitScore
                        
    
                blastResultFile_fh.close()
                
                bitScoreMatrixFile = NGS_Util.createFilePath(self.pfamBlastStaticticsDir, pfamID + "_AllvsALL_BlastBitScores.txt")
                
                self.writeMatrixToFile(bitScoreMatrix, bitScoreMatrixFile)
                
                return bitScoreMatrix
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def gettNormalizedBlastBitScoreMatrix(self, pfamID, bitScoreMatrix, blastResultFile):
    
        try:

            print "gettNormalizedBlastBitScoreMatrix"
            
                        
            normalizedBitScoreMatrix = self.createBlastBitScoreMatrix(blastResultFile)
            
            for sequenceId, matchingSequenceIDsList in bitScoreMatrix.iteritems():
        
                for matchingSequenceID in matchingSequenceIDsList:
                    
                    sum = bitScoreMatrix[sequenceId][sequenceId] + bitScoreMatrix[matchingSequenceID][matchingSequenceID] - bitScoreMatrix[sequenceId][matchingSequenceID]
                    
                    normalizedBitScoreMatrix[sequenceId][matchingSequenceID] = bitScoreMatrix[sequenceId][matchingSequenceID] / sum


            normalizedBitScoreMatrixFile = NGS_Util.createFilePath(self.pfamBlastStaticticsDir, pfamID + "_NormalizedBitScores.txt")
            self.writeMatrixToFile(normalizedBitScoreMatrix, normalizedBitScoreMatrixFile)
            
            return normalizedBitScoreMatrix
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def getNormalizeBlastBitScoresValuesSortedList(self, normalizeBlastBitScoresMatrix):
    
        try:

            print "getSortedNormalizeBlastBitScoresValuesList"
            
            matchingSequenceIDsDict     = {}
            normalizeBlastBitScoresList = []
            
            for sequenceId, matchingSequenceIDsDict in normalizeBlastBitScoresMatrix.iteritems():
        
                normalizeBlastBitScoresList = normalizeBlastBitScoresList + matchingSequenceIDsDict.values()
            
            return sorted(normalizeBlastBitScoresList)
        
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def getAbsoluteDeviationSortedList(self, median, valuesList):
    
        try:

            print "getAbsoluteDeviation"
                        
            for index in range(len(valuesList) ):
                
                valuesList[index] = abs(valuesList[index] - median)
        
            return sorted(valuesList)
            
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def getMedian(self, valuesList):
    
        try:

            print "getMedian"
                        
            length = len(valuesList)
            
            mid = length / 2
            
            remainder = length % 2
            
            if remainder == 1:
            
                return (valuesList[mid])
            
            else:
                
                return (valuesList[mid] + valuesList[mid+1] )/2
            
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def writeMatrixToFile(self, matrix, matrixFile):
    
        try:

            print "writeMatrixToFile"
            
            matrixFile_fh = open(matrixFile, "w")

            header = "SequenceID"

            sequenceIDsList = sorted(matrix.keys())
            
            
            for sequenceID in sequenceIDsList:
                
                header += "\t" + sequenceID
        
                
            matrixFile_fh.write(header+"\n")

            
            for sequenceId in matrix.iterkeys():
                
                line = sequenceId
        
                for matchingSequenceID in sequenceIDsList:
                    
                    line += "\t" + str(matrix[sequenceId][matchingSequenceID])
                    
                matrixFile_fh.write(line+"\n")


            matrixFile_fh.close()
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""



    def getPfamMedianAbsoluteDeviation(self, pfamID, blastResultFile):
    
        try:

            print "getPfamBlastStatictics"
            
            bitScoreMatrix = self.getBlastBitScoreMatrix(pfamID, blastResultFile)
            
            normalizedBitScoreMatrix = self.gettNormalizedBlastBitScoreMatrix(pfamID, bitScoreMatrix, blastResultFile)
            
            normalizedBitScoresList = self.getNormalizeBlastBitScoresValuesSortedList(normalizedBitScoreMatrix)
            
            normalizedBitScoresListMedian = self.getMedian(normalizedBitScoresList)
            
            
            normalizedBitScoresDeviationsList = self.getAbsoluteDeviationSortedList(normalizedBitScoresListMedian, normalizedBitScoresList)
            
            normalizedBitScoresListMAD = self.getMedian(normalizedBitScoresDeviationsList)

            return normalizedBitScoresListMedian, normalizedBitScoresListMAD
    
        except Exception:
            
            print traceback.print_exc()
            
        return ""

