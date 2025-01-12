#!/usr/bin/env python3
from operator import ge
import tarfile
from ERGO.database import db_session, engine
from ERGO.models.Genome import Genome
from ERGO.models.Feature import *
from csv import DictReader, DictWriter
from itertools import combinations
from tqdm import tqdm
import json
from multiprocessing import Pool, TimeoutError
import time

def bbh_id(orf, target_genome, percent_identity):
    s1 = db_session.query(Sim).filter(Sim.query_genome_id == orf.genome.id,
                                      Sim.subject_genome_id == target_genome.id,
                                      Sim.query_id == orf.id).order_by(Sim.percent_identity.desc()).first()
    if s1:
        s2 = db_session.query(Sim).filter(Sim.query_genome_id == target_genome.id,
                                          Sim.query_id == s1.subject_id,
                                          Sim.subject_genome_id == orf.genome.id).order_by(Sim.percent_identity.desc()).first()
        if s2:
            if s2.subject_id == s1.query_id and s2.percent_identity > percent_identity and s1.percent_identity > percent_identity:
                return s2.query

def feature_clusters_modified(source_genome: str, seed_list: List[str], genomes: List[str], 
                              p_score_threshold=1E-100, 
                              protein_coverage_cutoff=0.20, 
                              match_length_cutoff=0.20,
                              percent_identity=50,
                              ) -> dict:
    """
    Clusters features from genomes based upon the similarities in the similarity table with the provided
    cutoffs.
    :param genomes: List of the genomes to cluster features on
    :param evalue: E value cutoff, all considered sims will have an evalue less than this.
    :param protein_coverage_cutoff: The maximum amount of difference as a ratio allowed between two compared proteins.
    :return: Dictionary of {"clusters": [], "cluster_map": {}, "buckets": {}}
    """
    from ERGO.database import db_session
    clusters = []
    clustered = set()
    buckets = {}
    cluster_map = {}
    genome_metrics = {}

    source_genome = Genome.by(short_name=source_genome)
    seed_list = [ORF.by(name=o) for o in seed_list]
    genomes = [Genome.by(short_name=g) for g in genomes]

    # Create all combinations of genome.short_name s
    for cl in range(1, len([source_genome] + genomes) + 1):
        for bi in combinations([g.short_name for g in [source_genome] + genomes], cl):
            buckets[frozenset(bi)] = {"clusters": [], "gene_count": 0}

    for genome in [source_genome] + genomes:
        metrics = {}
        for k in buckets.keys():
            metrics[k] = 0
        genome_metrics[genome.short_name] = metrics

    for o in seed_list:
            if o.id in clustered:
                continue
            clustered.add(o.id)
            cluster = [o]
            orf_length = o.protein_length
            query_sig_domains = set([d.pssm_id for d in o.domains ])
            for g in genomes:
                for sim in o.all_sims(limit=10000, p_score_threshold=p_score_threshold, include_genome_ids=[g.id], display="internal"):
                    if sim.subject_id in clustered:
                        continue
                    if sim.percent_identity < percent_identity:
                        continue
                    if orf_length > sim.subject.protein_length:
                        match_length = (sim.query_stop - sim.query_start) + 1
                        dif = (orf_length - match_length) / orf_length
                    else:
                        match_length = (sim.subject_stop - sim.subject_start) + 1
                        dif = (sim.subject.protein_length - match_length) / sim.subject.protein_length
                    if dif > match_length_cutoff:
                        continue
                    overlaps = False
                    for c in cluster:
                        if sim.subject.overlaps(c):
                            overlaps = True
                            break
                    if overlaps:
                        continue
                    cluster.append(sim.subject)
                    clustered.add(sim.subject.id)


            clustered_genomes = frozenset([c.genome.short_name for c in cluster])
            [clustered.add(c.id) for c in cluster]
            for c in cluster:
                genome_metrics[c.genome.short_name][clustered_genomes] += 1
            cluster = [c.name for c in cluster]
            clusters.append(cluster)
            for c in cluster:
                cluster_map[c] = cluster
            buckets[clustered_genomes]["clusters"].append(cluster)
            buckets[clustered_genomes]["gene_count"] += len(cluster)
    
        
    db_session.remove()
    return {"clusters": clusters, "buckets": buckets, "cluster_map": cluster_map, "metrics": genome_metrics}


hspn = Genome.by(short_name="HSPN")
gene_list = {}

target_list = ["RTNK", "RJON"]
target_genomes = [Genome.by(short_name=g) for g in target_list]
genome_list = ["HSPN"] + target_list
cores = 40


genes = [f for f in db_session.query(Feature).filter_by(genome_id=hspn.id, feature_type="gene")]

jobs = []
with open("grouping.csv", "w") as outfh:
    writer = DictWriter(outfh, fieldnames=["star", "gene"] + genome_list)
    writer.writeheader()
    for gene in tqdm(genes, total=len(genes)):
        orfs = []
        for c in gene.children:
            if c.other_feature.feature_type == FeatureType.orf:
                orfs.append(c.other_feature.name)
        if len(orfs) == 0:
            print(F"No ORFs for {gene.name}.")
            continue
        jobs.append((gene.name, feature_clusters_modified, (genome_list[0],
                                                            orfs,
                                                            target_list,
                                                            .01,
                                                            0.40,
                                                            0.20,
                                                            60,
                                                            )))
    db_session.remove()
    db_session.close()
    engine.dispose()
    with Pool(processes=cores) as pool:
            running_jobs = []
            for (k, f, args) in jobs:
                running_jobs.append((k, pool.apply_async(f, args)))
            done = False
            pbar = tqdm(total=len(running_jobs))
            while len(running_jobs) > 0:
                remaning_jobs = []
                for (k, j) in running_jobs:
                    if j.ready():
                        res = j.get(timeout=1)
                        with open(F"{k}.group", "w") as fh:
                            metrics = dict(star="", gene=k)
                            for g in genome_list:
                                metrics[g] = []
                            for c in res["clusters"]:
                                for o in c:
                                    for g in genome_list:
                                        if g in o:
                                            metrics[g].append(o)
                            if len(metrics["RTNK"]) == 0:
                                metrics["star"] += "*"
                            if len(metrics["RJON"]) == 0:
                                metrics["star"] += "†"
                            for m in genome_list:
                                metrics[m] = ", ".join(metrics[m])
                            writer.writerow(metrics)
                            fh.write(json.dumps(res["clusters"], indent=4))
                        pbar.update(1)
                        pbar.set_description(k)
                    else:
                        remaning_jobs.append((k,j))
                running_jobs = remaning_jobs
                time.sleep(0.2)
            pbar.close()
