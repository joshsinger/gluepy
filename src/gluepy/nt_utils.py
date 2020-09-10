'''
Created on 8 Sep 2020

@author: joshsinger
'''

'''
          Nucleic Acid Code    Meaning                              Mnemonic
       ---------------------------------------------------------------------------------
          A                    A                                    Adenine
          C                    C                                    Cytosine
          G                    G                                    Guanine
          T                    T                                    Thymine
          U                    U                                    Uracil
          R                    A or G                               puRine
          Y                    C, T or U                            pYrimidines
          K                    G, T or U                            bases which are Ketones
          M                    A or C                               bases with aMino groups
          S                    C or G                               Strong interaction
          W                    A, T or U                            Weak interaction
          B                    not A (i.e. C, G, T or U)            B comes after A
          D                    not C (i.e. A, G, T or U)            D comes after C
          H                    not G (i.e., A, C, T or U)           H comes after G
          V                    neither T nor U (i.e. A, C or G)     V comes after U
          N                    A C G T U                            Nucleic acid
'''   
     
all_ambig_nts = "ACGTRYKMSWBDHVN"     
    
ambig_nt_to_concrete_nts = {}

ambig_nt_to_concrete_nts["A"] = "A"
ambig_nt_to_concrete_nts["C"] = "C"
ambig_nt_to_concrete_nts["G"] = "G"
ambig_nt_to_concrete_nts["T"] = "T"
ambig_nt_to_concrete_nts["R"] = "AG"
ambig_nt_to_concrete_nts["Y"] = "CT"
ambig_nt_to_concrete_nts["K"] = "GT"
ambig_nt_to_concrete_nts["M"] = "AC"
ambig_nt_to_concrete_nts["S"] = "CG"
ambig_nt_to_concrete_nts["W"] = "AT"
ambig_nt_to_concrete_nts["B"] = "CGT"
ambig_nt_to_concrete_nts["D"] = "AGT"
ambig_nt_to_concrete_nts["H"] = "ACT"
ambig_nt_to_concrete_nts["V"] = "ACG"
ambig_nt_to_concrete_nts["N"] = "ACGT"
