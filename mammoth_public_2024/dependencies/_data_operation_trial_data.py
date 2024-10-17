#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

import quantities as pq
import yaml
from SmartNeo.user_layer.dict_to_neo import templat_neo
import os

FILEPATH = os.path.dirname(os.path.abspath(__file__))

class _Path:
    pass

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
    
    def ObjectStatusRecordOperation(InputData,userElement):
        '''
        Transform userElement from .mat type to dict.
            The dict should be obey the rule from 'TransformationMap' in '.BehaviorBaseInterface'.
                The data under 'ObjectStatusRecord' is handled by 'irr' key word.
                    Because the time of the sample is given and the data is sampled irregularly without fixed sample rate

        Parameters
        ----------
        Same as 'Ml2InputData'

        Returns
        -------
        None.

        '''
        
        # initial data structure from the template
        InputData['ObjectStatusRecord']['Position']['irr'] = templat_neo['irr'].copy()
        InputData['ObjectStatusRecord']['Status']['irr'] = templat_neo['irr'].copy()
        
        # allocate appointed key
        # Under 'ObjectStatusRecord', there are two variable, 'Position' and 'Status'
        # all of them is saved with 'irr' key word
        
        # allocate value to 'signal'
        InputData['ObjectStatusRecord']['Position']['irr']['signal'] = [i[0].squeeze() for i in userElement['ObjectStatusRecord']['Position'][0][0]]*pq.dimensionless
        InputData['ObjectStatusRecord']['Status']['irr']['signal'] = list(userElement['ObjectStatusRecord']['Status'][0][0].squeeze())*pq.dimensionless
        
        # allocate value to 'times'
        InputData['ObjectStatusRecord']['Position']['irr']['times'] = list(userElement['ObjectStatusRecord']['Time'][0][0].squeeze())*pq.ms
        InputData['ObjectStatusRecord']['Status']['irr']['times'] = list(userElement['ObjectStatusRecord']['Time'][0][0].squeeze())*pq.ms
        
        # give a t_start
        InputData['ObjectStatusRecord']['Position']['irr']['t_start'] = 0*pq.ms
        InputData['ObjectStatusRecord']['Status']['irr']['t_start'] = 0*pq.ms
    
    
    def BehavioralCodesOperation(InputData,userElement):
        '''
        The data under 'BehavioralCodes' label contains time and event marker.
            This data type should be handled by 'event'

        Parameters
        ----------
        Same as 'Ml2InputData'

        Returns
        -------
        None.

        '''
        
        # initial data structure from the template
        InputData['BehavioralCodes']['event'] = templat_neo['event'].copy()
        
        # allocate data to appointed key
        InputData['BehavioralCodes']['event']['labels'] = list(userElement['BehavioralCodes']['CodeNumbers'][0][0].astype(float).squeeze())
        InputData['BehavioralCodes']['event']['times'] = list(userElement['BehavioralCodes']['CodeTimes'][0][0].squeeze())*pq.ms
    
    def OtherOperation(key):
        def WarpFunction(InputData,userElement):
            keyword = list(InputData[key].keys())[0]
            InputData[key][keyword] = userElement[key][0][0] if key!='AbsoluteTrialStartTime' else userElement[key][0][0]*pq.ms
        return WarpFunction
    
    def UserVarsOperation(InputData,userElement):
        for i in userElement['UserVars'].dtype.names:
            InputData['UserVars'][i] = {}
            InputData['UserVars'][i]['des'] = userElement['UserVars'][i][0][0]
            
    def AnalogDataOperation(InputData,userElement):
        ''' No data can be handled by "ana"'''
        
        pass
    
    # Build the map
    Template = yaml.safe_load(open(os.path.join(FILEPATH,'template_trial_data.yml')))
    dataOperationMap['ObjectStatusRecord'] = ObjectStatusRecordOperation
    dataOperationMap['BehavioralCodes'] = BehavioralCodesOperation
    dataOperationMap['UserVars'] = UserVarsOperation
    
    for i in set(Template.keys()).difference(['ObjectStatusRecord','BehavioralCodes','UserVars']):
        dataOperationMap[i] = OtherOperation(i)
    
    return dataOperationMap
        
        
        
        
        
        
