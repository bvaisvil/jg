#!/usr/bin/env python3

from multiprocessing import Pool
from subprocess import Popen, STDOUT
import argparse
import os
import tqdm

prank = "/mnt/big/aligners/prank-msa/src/prank"

def run_alignment(path) -> None:
    cmd = ["nice", prank, F"-d={path}/in.fasta", F"-o={path}/out", "-showall", "-codon", "-f=phylips"]
    print(cmd)
    cmd = Popen(cmd, stderr=STDOUT, stdout=open(f"{path}/prank.log", "w"), close_fds=True)
    result = cmd.wait()
    return (result, path)

def completed(result, pbar):
    if result[0] != 0:
        print(F"{result[1]} failed.")
    else:
        pbar.update(1)
        pbar.set_description(os.path.split(result[1])[-1])

def run(args):
    
    with Pool(args.processes) as pool:
        results = []
        pbar = tqdm.tqdm(desc="Running MSAS", total=len(os.listdir(args.dir)))
        for p in os.listdir(args.dir):
            print(p)
            results.append(pool.apply_async(run_alignment, args=(os.path.join(args.dir, p), ), callback=lambda r: completed(r, pbar)))
        for r in results:
            r.get()
        pbar.close()


parser = argparse.ArgumentParser(
                    prog = 'Run MSAS')
parser.add_argument("-d", "--dir", dest="dir", required=True)
parser.add_argument("-p", "--processes", dest="processes", required=False, default=20, type=int)
parser.set_defaults(func=run)
args = parser.parse_args()
args.func(args)