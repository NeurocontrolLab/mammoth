#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 16:12:28 2020

@author: cuilab
"""

import os
import yaml
from _data_operation import BuilddataOperationMap
from SmartNeo.user_layer.base_user_interface import UserInterface


FILEPATH = os.path.dirname(os.path.abspath(__file__))
# a class for .mat data transformation is built by inheriting the base class 'UserInterface'

class AIEShare(UserInterface):
    ''' The template data structure is saved to class variable '''

    _template = yaml.safe_load(open(os.path.join(FILEPATH,'template_bhv_data.yml')))
    _data_operation_map = BuilddataOperationMap()


class BRShare(UserInterface):
    ''' The template data structure is saved to class variable '''

    _template = yaml.safe_load(open(os.path.join(FILEPATH,'template_neural_data.yml')))
    _data_operation_map = BuilddataOperationMap()


class Ml2NeoTrial(UserInterface):
    ''' The template data structure is saved to class variable '''

    _template = yaml.safe_load(open(os.path.join(FILEPATH,'template_trial_data.yml')))
    _data_operation_map = BuilddataOperationMap()