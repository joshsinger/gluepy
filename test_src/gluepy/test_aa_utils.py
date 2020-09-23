'''
Created on 8 Sep 2020

@author: joshsinger
'''
import unittest
from gluepy.aa_utils import amb_tripl_to_result

class Test(unittest.TestCase):
    # Example 1
    # ambiguous NT triplet: YAY (Y = C/T)
    # concrete triplets set: CAC, CAT, TAC, TAT
    # AAs: Y (CAC, CAT), H (TAC, TAT)
    # Deleting CAC, CAT from the set gives TAC, TAT, equivalent to ambiguous triplet TAY, != YAY so Y is definite.
    # Deleting TAC, TAT from the set gives CAC, CAT, equivalent to ambiguous triplet CAY, != YAY so H is definite.
    def test1(self):
        result = amb_tripl_to_result['YAY']
        self.assertEquals(result.aa, 'X')
        self.assertEquals(result.concrete_triplets, ['CAC', 'CAT', 'TAC', 'TAT'])
        self.assertEquals(result.possible_aas, 'HY')
        self.assertEquals(result.definite_aas, 'HY')


    # Example 2
    # ambiguous NT triplet: TSR (S = C/G, R = A/G)
    # concrete triplets set: TCA, TCG, TGA, TGG
    # AAs: * (TGA), S (TCA, TCG), W (TGG)
    # Deleting TCA, TCG from the set gives TGA, TGG, equivalent to ambiguous triplet TGR, != TSR so S is definite.
    # Deleting TGA from the set gives TCA, TCG, TGG, equivalent to ambiguous triplet TSR so * is merely possible.
    # Deleting TGG from the set gives TCA, TCG, TGA, equivalent to ambiguous triplet TSR so W is merely possible.
    
    def test2(self):
        result = amb_tripl_to_result['TSR']
        self.assertEquals(result.aa, 'X')
        self.assertEquals(result.concrete_triplets, ['TCA', 'TCG', 'TGA', 'TGG'])
        self.assertEquals(result.possible_aas, '*SW')
        self.assertEquals(result.definite_aas, 'S')
        
    # Example 3
    # ambiguous NT triplet: ACN
    # concrete triplets: ACA, ACC, ACG, ACT
    # AAs: T (ACA, ACC, ACG, ACT)
    # deleting these is empty so T is definite.
    def test3(self):
        result = amb_tripl_to_result['ACN']
        self.assertEquals(result.aa, 'T')
        self.assertEquals(result.concrete_triplets, ['ACA', 'ACC', 'ACG', 'ACT'])
        self.assertEquals(result.possible_aas, 'T')
        self.assertEquals(result.definite_aas, 'T')
    
    # Example 4
    # ambiguous NT triplet: ANN
    # concrete triplets: AAA, AAC, AAG, AAT,   ACA, ACC, ACG, ACT,  AGA, AGC, AGG, AGT,   ATA, ATC, ATG, ATT
    # AAs: KNTRSIM
    # Looking at T (ACA, ACC, ACG, ACT), deleting these leaves a set equivalent to ADN (D=A/G/T),
    # This is different from ANN at position 2 however there is an N in the original so T is not definite.
    # T was "relying" on the N to provide its C at position 2.
    def test4(self):
        result = amb_tripl_to_result['ANN']
        self.assertEquals(result.aa, 'X')
        self.assertEquals(result.concrete_triplets, ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT'])
        self.assertEquals(result.possible_aas, 'IKMNRST')
        self.assertEquals(result.definite_aas, '')




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()