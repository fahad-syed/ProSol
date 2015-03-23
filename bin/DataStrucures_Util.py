#!/usr/bin/env python

import os
import sys
import subprocess
import operator
import traceback
import time


def convertDictionaryListValuesToSetValues(dictionary):
    
    try:

        for key, value in dictionary.iteritems():
            
            dictionary[key] = set(value)
        
        return dictionary
      
    except Exception:
        
        print traceback.print_exc()
        
    return None


def addUniqueListValueToDictionary(key, value, dictionary):
    
    try:

        if not dictionary.has_key(key):
            
            dictionary[key]=[value]
            
        else:
            if dictionary[key].count(value) < 1:
                dictionary[key].append(value)

        return dictionary
      
    except Exception:
        
        print traceback.print_exc()
        
    return None


def addListValueToDictionary(key, value, dictionary):
    
    try:

        if not dictionary.has_key(key):
            
            dictionary[key]=[value]
            
        else:
            
            dictionary[key].append(value)

        return dictionary
      
    except Exception:
        
        print traceback.print_exc()
        
    return None


def convertDictionaryListValuesToSetValues(dictionary):
    
    try:

        for key, value in dictionary.iteritems():
            
            dictionary[key] = set(value)
        
        return dictionary
      
    except Exception:
        
        print traceback.print_exc()
        
    return None


def getDictionaryKeysAsTabSeperatedValues(dictionary):
    
    try:
        
        line = ""

        for key in dictionary.iterkeys():
            
            line += key  + "\t"
        
        return line[0:len(line)-1]
      
    except Exception:
        
        print traceback.print_exc()
        
    return None


def writeDictionaryListValuesToFile(dictionary,file):
    
    try:

        dict_fh = open(file, "w")
        
        line = ""
        vals = ""
        
        for key, values in dictionary.iteritems():

            vals = ""
            
            for val in values:
                
                vals += str(val)  + "\t"
                            
            if len(vals) > 0:
                
                line = key + "\t" + vals[0:len(vals)-1] + "\n"
                
            else:
                
                line = key + "\n"
        
            dict_fh.write(line)

        dict_fh.close()
      
    except Exception:
        
        print traceback.print_exc()
        

def writeDictionaryListValuesToFileWithHeader(header, dictionary,file):
    
    try:

        dict_fh = open(file, "w")
        
        line = ""
        vals = ""
        
        dict_fh .write(header+"\n")
        
        for key, values in dictionary.iteritems():
            
            vals = ""
            
            for val in values:
                
                vals += str(val) + "\t" 
            
            if len(vals) > 0:
                
                line = key + "\t" + vals[0:len(vals)-1] + "\n"
                
            else:
                
                line = key + "\n"
        
            dict_fh.write(line)

        dict_fh.close()
      
    except Exception:
        
        print traceback.print_exc()
        
        


def mergeDictionaries(dict1, dict2):
    
    try:
        
        dictionary = {}

        for key, value in dict1.iteritems():
            
            dictionary[key] = [value]
            
        for key, value in dict2.iteritems():
            
            if dictionary.has_key(key):
                
                dictionary[key].append(value)
                
            else:
                
                dictionary[key] = [value]
        
        return dictionary
      
    except Exception:
        
        print traceback.print_exc()
        
    return None
        
        

def mergeListsToUniqueValuesList(list1, list2):
    
    try:

        resulting_list = sorted ( list ( set(list1 + list2) ) )
        
        return resulting_list
      
    except Exception:
        
        print traceback.print_exc()
        
    return None



def writeListValuesToFileWithHeader(header, list,file):
    
    try:

        list_fh = open(file, "w")
        
        line = ""
        
        list_fh .write(header+"\n")
        
        for val in list:

            line += str(val) + "\t" 
                
        line += "\n"
        
        list_fh.write(line)

        list_fh.close()
      
    except Exception:
        
        print traceback.print_exc()
        
        