#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 14:04:07 2022

@author: lfleming
"""

import unittest
import calc

class TestCalc(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        print('SetUpClass\n')
        
    @classmethod
    def tearDownClass(self):
        print('TearDownClass\n')
        
        
    def setUp(self):
        print('setUp')
        
    def tearDown(self):
        print('tearDown\n')
    
    def test_add(self):
        print('test_add')
        self.assertEqual(calc.add(5, 10), 15)
        self.assertEqual(calc.add(-5, 10), 5)
        self.assertEqual(calc.add(0, 10), 10)
        self.assertEqual(calc.add(-5, -10), -15)
        
        
    def test_subtract(self):
        print('test_subtract')
        self.assertEqual(calc.subtract(5, 10), -5)
        self.assertEqual(calc.subtract(-5, 10), -15)
        self.assertEqual(calc.subtract(0, 10), -10)
        self.assertEqual(calc.subtract(-5, -10), 5)
        
        
    def test_multiply(self):
        print('test_multiply')
        self.assertEqual(calc.multiply(5, 10), 50)
        self.assertEqual(calc.multiply(-5, 10), -50)
        self.assertEqual(calc.multiply(0, 10), 0)
        self.assertEqual(calc.multiply(-5, -10), 50)
        
    def test_divide(self):
        print('test_divide')
        self.assertEqual(calc.divide(2, 1), 2)
        self.assertEqual(calc.divide(2, -1), -2)
        self.assertEqual(calc.divide(-2, -1), 2)
        self.assertEqual(calc.divide(-2, 1), -2)
        self.assertEqual(calc.divide(5, 2), 2.5)
        self.assertRaises(ValueError, calc.divide, 10, 0)
        
        with self.assertRaises(ValueError):
            calc.divide(10, 0)
                
                
if __name__ == '__main__':
    unittest.main()