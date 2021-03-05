#!/usr/bin/env python3

import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/../')
import unittest
from unittest.mock import patch
from unittest import mock
import filecmp
import numpy as np

import pmf

run_dir = os.path.dirname(os.path.realpath(__file__))+'/'
class TestSum(unittest.TestCase):
    def test_set_to_zero(self):
        e1 = np.array([20, 10])
        e2 = np.array([-20, -10])
        result = pmf.set_to_zero(e1)
        print(result)
        self.assertEqual(list(result), [10,0])
        result = pmf.set_to_zero(e2)
        self.assertEqual(list(result), [-10,0] )



if __name__ == '__main__':
    unittest.main()
