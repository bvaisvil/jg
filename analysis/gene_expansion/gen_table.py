from ERGO.database import db_session
from ERGO.models.Feature import Feature
import csv

with open("grouping.csv") as fh:
    line_number = 0
    header = []
    for line in csv.DictReader(fh):    
        f = Feature.by(name=line["gene"])
        counts = {}
        orfs = {}
        for k in line.keys():
            if k in ["gene", "star"]:
                continue
            if line_number == 0:
                header.append(k)
            orfs[k] = line[k].split(", ")
            
            overlap_sets = set()
            for o1 in orfs[k]:
                overlap_set = set([o1])
                for o2 in orfs[k]:
                    if o1 == o2:
                        continue
                    f1 = Feature.by(name=o1)
                    f2 = Feature.by(name=o2)
                    if f1.overlaps(f2):
                        overlap_set.add(o2)
                overlap_sets.add(frozenset(overlap_set))
            counts[k] = len(overlap_sets)
                        
        if line_number == 0:
            print("gene\t" + "\t".join(header))
        print(F"{f.gene_name}\t" + "\t".join([str(counts[h]) for h in header]))
        
        line_number += 1