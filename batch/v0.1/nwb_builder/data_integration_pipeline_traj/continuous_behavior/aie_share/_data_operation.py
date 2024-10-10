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
    
    @NullLabel('AnalogData')
    def AnalogDataOperation(InputData,userElement):
        InputData['AnalogData'] = userElement['AnalogData']
        
    @NullLabel('IrregularSampledData')
    def IrregularSampledDataOperation(InputData,userElement):
        InputData['IrregularSampledData'] = userElement['IrregularSampledData']
    
    @NullLabel('Event')
    def EventOperation(InputData,userElement):
        InputData['Event'] = userElement['Event']
    
    @NullLabel('name')
    def nameOperation(InputData,userElement):
        InputData['name'] = userElement['name']
    
    # Build the map
    dataOperationMap['AnalogData'] = AnalogDataOperation
    dataOperationMap['IrregularSampledData'] = IrregularSampledDataOperation
    dataOperationMap['Event'] = EventOperation
    dataOperationMap['name'] = nameOperation
    
    return dataOperationMap
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
