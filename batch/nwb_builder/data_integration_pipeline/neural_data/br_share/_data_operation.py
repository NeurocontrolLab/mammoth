#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""


def NullLabel(label:str):
    def CheckNull(func):
        def warpped_function(InputData,userElement):
            if isinstance(userElement[label],str):
                if userElement[label]=='null':
                    InputData[label] = 'null'
                    return
            func(InputData,userElement)
        return warpped_function
    return CheckNull

def BuilddataOperationMap()->dict:
    '''
    A map to handle .mat data given by monkeylogic

    Returns
    -------
    dict
        A dict contains column labels and function for operation.

    '''
    
    # initialize the map (dict)
    dataOperationMap = {}
    
    @NullLabel('RecordingSystemEvent')
    def RecordingSystemEventOperation(InputData,userElement):
        InputData['RecordingSystemEvent']['event'] = userElement['RecordingSystemEvent']
        
    @NullLabel('LFP')
    def LFPOperation(InputData,userElement):
        InputData['LFP']['irr'] = userElement['LFP']
    
    @NullLabel('SBP')
    def SBPOperation(InputData,userElement):
        InputData['SBP']['irr'] = userElement['SBP']
    
    @NullLabel('Spike')
    def SpikeOperation(InputData,userElement):
        InputData['Spike'] = userElement['Spike']
    
    @NullLabel('TCR')
    def TCROperation(InputData,userElement):
        InputData['TCR'] = userElement['TCR']
    
    @NullLabel('name')
    def nameOperation(InputData,userElement):
        InputData['name'] = userElement['name']
    
    # Build the map
    dataOperationMap['RecordingSystemEvent'] = RecordingSystemEventOperation
    dataOperationMap['LFP'] = LFPOperation
    dataOperationMap['Spike'] = SpikeOperation
    dataOperationMap['SBP'] = SBPOperation
    dataOperationMap['TCR'] = TCROperation
    dataOperationMap['name'] = nameOperation
    
    return dataOperationMap
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
