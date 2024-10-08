#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 22:03:01 2021

@author: cuilab
"""

from SmartNeo.user_layer.base_user_interface import UserInterface
import yaml
from _data_operation import BuilddataOperationMap
import os

FILEPATH = os.path.dirname(os.path.abspath(__file__))
# a class for .mat data transformation is built by inheriting the base class 'UserInterface'
class Ml2NeoTrial(UserInterface):
    ''' The template data structure is saved to class variable '''

    _template = yaml.safe_load(open(os.path.join(FILEPATH,'template.yml')))
    _data_operation_map = BuilddataOperationMap()
        
        
        
        
        
        
        
