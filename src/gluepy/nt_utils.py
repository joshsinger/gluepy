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
     
ambig_nt_to_concrete_nts = {}

ambig_nt_to_concrete_nts["A"] = ["A"]
ambig_nt_to_concrete_nts["C"] = ["C"]
ambig_nt_to_concrete_nts["G"] = ["G"]
ambig_nt_to_concrete_nts["T"] = ["T"]
ambig_nt_to_concrete_nts["R"] = ["A", "G"]
ambig_nt_to_concrete_nts["Y"] = ["C", "T"]
ambig_nt_to_concrete_nts["K"] = ["G", "T"]
ambig_nt_to_concrete_nts["M"] = ["A", "C"]
ambig_nt_to_concrete_nts["S"] = ["C", "G"]
ambig_nt_to_concrete_nts["W"] = ["A", "T"]
ambig_nt_to_concrete_nts["B"] = ["C", "G", "T"]
ambig_nt_to_concrete_nts["D"] = ["A", "G", "T"]
ambig_nt_to_concrete_nts["H"] = ["A", "C", "T"]
ambig_nt_to_concrete_nts["V"] = ["A", "C", "G"]
ambig_nt_to_concrete_nts["N"] = ["A", "C", "G", "T"]
