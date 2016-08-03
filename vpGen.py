# -*- coding: utf-8 -*-
"""
Created on Sun May 15 2016

@author: Matthew Carse
"""
#@ The vpGen class takes an anchor and a chain residue and creates an identifier from these residues.
#@ When the featureGen function is called with an appropriate peptide sequence, the virtual potential for the
#@ sequence is calculated using the anchor and chain residues. The virtual potential is displayed as a 6-item 
#@ feature set, with each 'item' corresponding to a range of 0.5 from 0 to 3. The maximum value of a virtual
#@ potential is 2.92896825397.

class vpGen:     
    # constructor taking anchor and chain residues and creating identifier
    def __init__(self, anchor, chain):
        # use join for instance variables as anchor/chain are lists due to .split() operation
        self.anchor = ''.join(anchor)
        self.chain = ''.join(chain)
        self.id = self.anchor+self.chain  
    
    def featureGen(self, seq):
        "Derive virtual potentials given a sequence, anchor residue and chain residue" 
        # find the position of the first anchor residue in the sequence (has to be at least index 10 due to window size of 10)
        ancPos = seq.find(self.anchor,10)
        
        if ancPos == -1:
            # default vp set
            # float values for future calculations
            vps = [-1.,-1.,-1.,-1.,-1.,-1.]
            return vps
        else:
            # recursively call featureGen until anchor not found in remaining sequence
            vps = self.featureGen(seq[ancPos-9:])
            
            # if first iteration with anchor existing, set virtual potential counts to 0
            # float values for future calculations
            if vps == [-1.,-1,-1.,-1.,-1.,-1.]:
                vps = [0.,0.,0.,0.,0.,0.]
                
            # reverse sequence for ease of calculation using indices
            window = seq[(ancPos-10):ancPos][::-1]
            
            # iterate through the 10-residue window, searching for chain residue
            vp = 0
            for i, c in enumerate(window):
                if c==self.chain: 
                    # if chain residue is found, calculate virtual potential
                    vp += 1.0/(i+1)  

            # look up correct index given virtual potential and increment 'histogram'
            if vp > 0:
                vps[self.vpDict(vp)]+=1

            return vps
        
    def vpDict(self, vp):
        "Returns array index corresponding to virtual potential value range"
        if (vp >= 0.0) & (vp < 0.5): return 0
        if (vp >= 0.5) & (vp < 1.0): return 1
        if (vp >= 1.0) & (vp < 1.5): return 2
        if (vp >= 1.5) & (vp < 2.0): return 3
        if (vp >= 2.0) & (vp < 2.5): return 4
        if (vp >= 2.5) & (vp < 3.0): return 5