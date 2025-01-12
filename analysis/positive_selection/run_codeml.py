#!/usr/bin/env python3

from multiprocessing import Pool
from subprocess import Popen, STDOUT
import argparse
import os
import tqdm
from typing import Tuple, List
from Bio.Phylo.PAML import codeml
import re
import shutil

codeml_bin = "/mnt/big/tools/paml-4.10.6/src/codeml"
reorder_alignments = "reorder_alignments.py"
rscript = "/home/benjamin/miniconda3/bin/Rscript"
convert_tree_R = "convert_tree.R"
tort_re = re.compile(r"(TORT:(?:\d|\.)+)")

def run_codeml(template, path) -> Tuple[int, str]:
    shutil.copy(template, os.path.join(path, "codeml.ctl"))
    cwd = os.getcwd()
    os.chdir(path)
    cmd = Popen(("nice", codeml_bin), stderr=STDOUT, stdout=open(f"{template}.log", "w"), close_fds=True)
    result = cmd.wait()
    if result == 0 and os.path.exists(template + ".result"):
        parsed = codeml.read(template + ".result")
    else:
        parsed = None
    os.chdir(cwd)
    return (result, path, template, parsed)

def convert_tree(path):
    cwd = os.getcwd()
    os.chdir(path)
    cmd = Popen(("nice", rscript, convert_tree_R))
    result = cmd.wait()
    os.chdir(cwd)
    if result != 0 or not os.path.exists(os.path.join(path, "out.best.unrooted.dnd")):
        raise RuntimeError("Couldn't convert tree for: {path}")
        

def mark_tort(match_obj):
    if match_obj.group(1) is not None:
        return match_obj.group(1) + " #1"

def run_codemls(path, controls) -> List[Tuple[int, str]]:
    results = []
    if not os.path.exists(os.path.join(path, "out.best.dnd")):
        print(F"No tree for {path}")
        return []

    if "codeml.m2a" in controls:
        try:
            convert_tree(path)
        except RuntimeError as e:
            print(F"Skipping {path}...")
            print(e)
            return results
        with open(os.path.join(path, "out.best.unrooted.dnd")) as fh:
            with open(os.path.join(path, "out.best.marked.dnd"), "w") as outfh:
                for line in fh:
                    if "TORT" in line:
                        line = re.sub(tort_re, mark_tort, line)
                    outfh.write(line)
    if "codeml.m8" in controls:
        cmd = Popen(("python3", reorder_alignments, "-i", F"{path}/out.best.phy", "-o", f"{path}/out.best.reordered.phy", "-f", "HUMAN"))
        cmd.wait()
    for t in controls:
        (result, p, t, rob) = run_codeml(t, path)
        if result == 0 and rob is not None:
            if "NSsites" not in rob:
                raise ValueError(rob)
            if t == "codeml.m1":
                if 1 not in rob["NSsites"]:
                    raise ValueError(rob["NSsites"])
                if "lnL" not in rob["NSsites"][1]:
                    raise ValueError(rob["NSsites"][1])
                if "parameters" not in rob["NSsites"][1]:
                    raise ValueError(rob["NSsites"][1])
                results.append([p, 
                rob['NSsites'][1]['lnL'], 
                rob['NSsites'][1]['parameters']['site classes'][0]['omega'],
                rob['NSsites'][1]['parameters']['site classes'][1]['omega'],
                ])
            elif t == "codeml.m2a":
                results.append([p,
                                rob['NSsites'][2]['lnL'],
                                rob['NSsites'][2]['parameters']['site classes'][2]['branch types']['background'],
                                rob['NSsites'][2]['parameters']['site classes'][2]['branch types']['foreground'],
                                rob['NSsites'][2]['parameters']['site classes'][3]['branch types']['background'],
                                rob['NSsites'][2]['parameters']['site classes'][3]['branch types']['foreground'],
                                ])
            elif t == "codeml.m7":
                results.append([p,
                                rob['NSsites'][7]['lnL'],
                                ])
                print(results)
            elif t == "codeml.m8":
                results.append([p,
                                rob['NSsites'][8]['lnL'],
                                rob['NSsites'][8]['parameters']['w']
                                ])
                print(results)
    return results

def completed(results, path, pbar):
    pbar.update(1)
    if len(results) > 0 and len(results[0]) >0:
        print(results)
        pbar.set_description(results[0][0])
    else:
        pbar.set_description(path)

def run(args):
    analyze_these = set()
    filter = False
    if args.filter:
        filter = True
        with open(args.filter, "r") as fh:
            for line in fh:
                analyze_these.add(line.rstrip())
    with open(args.output, "w") as outfh:
        with Pool(args.processes) as pool:
            results = []
            pbar = tqdm.tqdm(desc="Running CodeMLs", total=len(os.listdir(args.dir)))
            for p in os.listdir(args.dir):
                parts = p.split("_")
                if filter and parts[1] not in analyze_these:
                    continue
                results.append(pool.apply_async(run_codemls, args=(os.path.join(args.dir, p), args.controls), callback=lambda r: completed(r, p, pbar)))
            for r in results:
                result = r.get()
                outfh.write(F"{p}\t")
                for line in result:
                    outfh.write("\t".join([str(f) for f in line]) + "\t")
                outfh.write("\n")
            pbar.close()


parser = argparse.ArgumentParser(
                    prog = 'Run CodeML')
parser.add_argument("-f", "--filter", dest="filter", required=False, default=None)
parser.add_argument("-d", "--dir", dest="dir", required=True)
parser.add_argument("-p", "--processes", dest="processes", required=False, default=20, type=int)
parser.add_argument("-o", "--output", dest="output", required=True)
parser.add_argument("-c", "--controls", dest="controls", default=["codeml.m1", "codeml.m2a"], nargs="+")
parser.set_defaults(func=run)
args = parser.parse_args()
args.func(args)