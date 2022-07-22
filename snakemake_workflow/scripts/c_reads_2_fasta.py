import pandas as pd
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

outfile = open(outfile, "w")
for chonk in pd.read_table(infile, chunksize = 1000):
    

    for index, row in chonk.iterrows():
        outfile.write(">" + "sample_id:" + 
                      row["sample_id"] + "_" + 
                      "totalumis:" + str(row["total_umis"]) + "\n" + row["c_sequence"] + "\n")

outfile.close()
