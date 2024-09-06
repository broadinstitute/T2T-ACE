#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 14:01:32 2022

@author: lfleming
"""

def add(x,y):
     return x+y
 
def subtract(x,y):
    return x-y

def multiply(x,y):
    return x*y

def divide(x,y):
    
    if y == 0:
        raise ValueError('Cannot divide by zero!')
    return x / y

