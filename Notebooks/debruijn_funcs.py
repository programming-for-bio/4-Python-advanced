#!/usr/bin/env python

"""
Functions for assembling contigs from sequenced reads using kmers
and a debruijn graph.
"""


import random
#from eulerian import eulerian_path


def random_sequence(nbases):
    return "".join((random.choice("ACGT") for i in range(nbases)))



def get_kmers(target, k):
    "returns k-mers dict for a string target"
    kmers = {}
    for i in range(0, len(target) - k + 1):
        kmer = target[i:i+k]
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
    return kmers



def get_reads(string, nreads, rlen):
    "returns nreads of len rlen drawn from string"  
    last_start = len(string) - rlen
    startpoints = [random.randint(0, last_start) for i in range(nreads)]
    reads = [string[i:i+rlen] for i in startpoints]
    return reads



def reads_to_kmers(reads, k):
    "stores kmers to dict uses update to join together kmer dict keys"
    kmers = {}
    for read in reads:
        ikmers = get_kmers(read, k)
        kmers.update(ikmers)
    return kmers
    
    

def get_debruijn_edges(kmers):
    " doc string"
    edges = set()
    kmers = tuple(kmers.keys())
    for k1 in kmers:
        for k2 in kmers:
            ## if xaa = aax then add (aa, ax)
            if k1[1:] == k2[:-1]:
                edges.add((k1[1:], k2[1:]))
    return edges



if __name__ == "__main__":
    
    # run a test
    random.seed(123)
    target = random_sequence(100)
    reads = get_reads(target, 200, 20)
    kmers = reads_to_kmers(reads, 8)
    edges = get_debruijn_edges(kmers)
    
    print(edges)
    