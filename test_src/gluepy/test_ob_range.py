'''
Created on 11 Sep 2020

@author: joshsinger
'''
import unittest

from gluepy.ob_range import OneBasedRange, OneBasedRangeException


class Test(unittest.TestCase):


    def testObRange1(self):
        self.assertEquals(OneBasedRange(1, 4).pull_from_str("ACTGGTA"), "ACTG")

    def testObRange2(self):
        with self.assertRaises(OneBasedRangeException):
            OneBasedRange(1, 10).pull_from_str("ACTGGTA")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()