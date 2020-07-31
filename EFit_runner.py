import ROOT
from rat import dsreader
import sys
import EFit

#outDir = "/SNO+/pdfs/" # for local VM
outDir = "/home/eigenboi/pdfs/" # for cedar
outFileName = "example.pdf"
index = 1

if __name__ == '__main__' :
    
    EFit.BiPoEnergy()