# Alexander M. Procton / 2.18.2018
# Based on De Bruijn functions by Deren Eaton for Programming and Data Science for Biologists

class Assembler5():
    """
    An object for constructing a debuijn graph from kmers of random reads
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
        self.target = "".join((random.choice("ACGT") for i in range(target_length)))
        
        
    def _get_reads(self, nreads, rlen):
        "returns nreads of len rlen drawn from string"  
        last_start = len(self.target) - rlen
        startpoints = [random.randint(0, last_start) for i in range(nreads)]
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
        "stores kmers to dict uses update to join together kmer dict keys"
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
          
        
    ## public functions
    def run(self, nreads, rlen, k):
        "generates reads and breaks them into kmers"
        self._get_reads(nreads, rlen)
        self._reads_to_kmers(k)
        self._get_debruijn_edges()
        self._get_eulerian_path()
        
    def test(self, target_size = [200, 500, 1000], nreads = [500, 1000, 5000], k = [10, 20, 30]):
        """
        Tests assembler with different values of target_size, nreads, and k
        Default values: target_size = [200, 500, 1000], nreads = [500, 1000, 5000], k = [10, 20, 30]
        """
        result_dict = {}
        
        for size in target_size:
            data = Assembler4(size)
            
            for n in nreads:
                for klen in k:
                    data.run(n, rlen=50, k=klen)

                    ## store result in dict
                    result = data.assembly == data.target
                    result_dict[(size, n, klen)] = result
                    
        return result_dict
    
    def plot(self):
        "Plots the De Bruijn graph after the edges have been obtained"
        e0 = [i[0] for i in self.edges]
        e1 = [i[1] for i in self.edges]
        toyplot.graph(e0, e1, tmarker=">", vlstyle={'font-size': '8px'});