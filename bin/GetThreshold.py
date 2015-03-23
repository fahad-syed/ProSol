#! /usr/bin/python

from sys import stdin,argv
from getopt import getopt
from math import *


def getAvgDegree(BlastFile,thresh):


	#degree
	degreeDict = {}

	#edges
	edgeDict = {}


	for input in file(BlastFile):
    		line = input.strip()
		source,target,prob = line.split(' ')
	
	
		#check for duplicate edges
        	if target in edgeDict.keys():
            		if source in edgeDict[target]:
                		continue

		if source in edgeDict.keys():
        		edgeDict[source].append(target)
	
		else:
			edgeDict[source] = [target]	

		#calculate edge-weight
        	if prob[0] == 'e':
            		prob = '1'+prob

        	prob = eval(prob)

        	if prob == 0.0:
            		continue

        	prob = int((-log10(prob)))

		if(prob < thresh):
			continue

		if source in degreeDict.keys():
			degreeDict[source].append(prob)
		else:
			degreeDict[source] = [prob]

	
		if target in degreeDict.keys():
			degreeDict[target].append(prob)
		else:
			degreeDict[target] = [prob]

	degree_sum = 0.0

	for weightList in degreeDict.values():

		degree = float(sum(weightList))
		degree_sum += degree


	return degree_sum / float(len(degreeDict.values()))






def get_threshold(BlastFile):

	change_Nsv = 0.0
	old_Nsv = getAvgDegree(BlastFile,0)

	
	
	

	for i in range(1,200):

		Nsv = getAvgDegree(BlastFile,i)

		
		change_Nsv = Nsv - old_Nsv

	
		if change_Nsv > 0:
			return "Threshold Found: %d"%i

		old_Nsv = Nsv
	
	return "No Threshold Found"

print get_threshold(argv[1])

