#!/usr/bin/env python


import toyplot
import random
import debruijn_funcs
from eulerian import eulerian_path

class Assembler6():
    """
    An object for constructing a debuijn graph from kmers of random 
reads
    """
    def __init__(self, target_length, random_seed=123):

        ## store attributes
        self.target = None 
        self.reads = None
        self.kmers = None
        self.edges = None
        self.assembly = None
        
        ## run init functions
        random.seed(random_seed)
        self._random_sequence(target_length)
        
        
    ## private functions
    def _random_sequence(self, target_length):
        self.target = "".join((random.choice("ACGT") for i in 
range(target_length)))
        
        
    def _get_reads(self, nreads, rlen):
        "returns nreads of len rlen drawn from string"  
        last_start = len(self.target) - rlen
        startpoints = [random.randint(0, last_start) for i in 
range(nreads)]
        self.reads = [self.target[i:i+rlen] for i in startpoints]
        
        
    def _get_kmers(self, string, k):
        "returns k-mers dict for a string target"
        kmers = {}
        for i in range(0, len(string) - k + 1):
            kmer = string[i:i+k]
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1
        return kmers
    

    def _reads_to_kmers(self, k):
        "stores kmers to dict uses update to join together kmer dict 
keys"
        kmers = {}
        for read in self.reads:
            ikmers = self._get_kmers(read, k)
            kmers.update(ikmers)
        self.kmers = kmers
        


    def _get_debruijn_edges(self):
        " return edges of the debruijn graph for a set of kmers"
        edges = set()
        kmers = tuple(self.kmers.keys())
        for k1 in kmers:
            for k2 in kmers:
                ## if xaa = aax then add (xaa, aax)
                if k1[1:] == k2[:-1]:
                    edges.add((k1, k2))
        self.edges = edges

        
        
    def _get_eulerian_path(self):
        """
        returns eulerian path through kmers joined as a string.
        Uses the loaded 'eulerian_path()' function from eulerian.py
        """
        try:
            epath = eulerian_path(self.edges)
            path = epath[0]
            for kmer in epath[1:]:
                path += kmer[-1]
            self.assembly = path
        except Exception:
            self.assembly = ""
          
        
    def _plot(self):
        e0 = [i[0] for i in self.edges]
        e1 = [i[1] for i in self.edges]
        toyplot.graph(e0, e1, tmarker=">", vlstyle={'font-size': 
'8px'});
    
    
    ## public function
    def run(self, nreads, rlen, k):
        "generates reads, breaks them into kmers, and generates a plot"
        self._get_reads(nreads, rlen)
        self._reads_to_kmers(k)
        self._get_debruijn_edges()
        self._get_eulerian_path()
        self._plot()              #automatically generate plot when 
running obj.run()
        
        
    def test(self, nreads = [500, 1000, 5000], kmer = [10, 20, 30]): #i 
want to allow the user to select their own test parameters, hopefully 
this works
        "test different parameters for effect on assembly"
        result_dict = {}
                   
        
        for i in [500, 1000, 5000]:     #for each input number of reads,
            for j in [10, 20, 30]:   #analyze each read by input kmer 
parameter
                self.run(nreads = i, rlen = 50, k = j)
                
                #store result
                result = self.assembly == self.target
                result_dict[(i,j)] = result    # i = number of reads, j 
= kmer parameter
                
            return result_dict
        
    def plot(self, nreads, rlen, k): #can't figure out how to make this 
work but might be cool
        self._plot() #nonetype objects are not iterable? need to put 
silent fxns inside this? maybe fix later
