import sys
import os
import pandas



def createopenbed(infile, outfile):
    df = pandas.read_csv(infile, sep="\t", names=["chr", "start", "stop", "name", "u1", "u2", "openstart", "openstop", "color", "n_subregions", "len_subregions", "startssubregions", "score", "u3", "u4"])
    newdf = df[["chr", "openstart", "openstop", "name", "score", "u1"]]
    newdf = newdf[newdf["openstart"]!=0]
    newdf.to_csv(outfile, sep="\t", index=False, header=False)



if __name__=="__main__":
        infile = sys.argv[1]
        outfile = sys.argv[2]
        print ("collecting open peaks ", infile, outfile)
        createopenbed(infile, outfile)
~                                         
