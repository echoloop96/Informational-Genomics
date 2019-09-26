import Bio
import Bio.Seq
from Bio import SeqIO
import Bio.Alphabet
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from statistics import mean, stdev
import math
import random
import statistics
import matplotlib.pyplot as plt
from pysuffixarray.core import SuffixArray
def nuc_alphabet():
    return ["A","C","G","T"]
'''

|----------------------------------------------------------------------------------|
| W: This framework is case sensitive please input your sequences all in uppercase |
|----------------------------------------------------------------------------------|


'''

def rand_seq(k):
    ''' Returns a random sequence with length k of generic_dna
        -----------
        Parameters:
            K(int) : sequence length
        --------
        Returns:
            Seq(Bio.Seq) : random sequence of length k
        
    
    '''
    seq="";
    
    for i in range (0,k) :
        r=random.randint(0,4)
        if r==0:
            seq=seq+"A"
        if r==1:
            seq=seq+"C"
        if r==2:
            seq=seq+"G"
        if r==3:
            seq=seq+"T"
    ret=Seq(seq,generic_dna)
    return ret;

def rand_seq_rand_length(k):
    
    ''' Returns a random sequence with a random length between 1 and k of generic_dna
        -----------
        Parameters:
            k(int) : length max of sequence
        --------
        Returns:
            Seq(Bio.Seq) : random Seq of length random between 1 and k
    '''
    
    seq="";
    length=random.randint(1,k)
    for i in range (0,length) :
        r=random.randint(0,4)
        if r==0:
            seq=seq+"A"
        if r==1:
            seq=seq+"C"
        if r==2:
            seq=seq+"G"
        if r==3:
            seq=seq+"T"
    ret=Seq(seq,generic_dna)
    return ret;

def get_alphabet(Seq):
    
    ''' Returns the alphabet of the sequence
        -----------
        Parameters:
            Seq(Bio.Seq) : reference sequence
        Returns:
            set(str) : set of strings (alphabet)
    '''
    
    alphabet = set()
    for i in Seq:
        alphabet.add(str(i))
    return alphabet

def print_kmers(Seq, k):
    
    ''' Prints the set of k-mers that occur in a given sequence 
        ----------
        Paramaters:
            Seq(Bio.Seq)
            K(int)  
    '''
    for i in range(len(Seq) - k +1):
        print(Seq[i:i+k])
    
def get_kmers(Seq, k):
    
    ''' Returns the set of k-mers that occur in a given sequence 
        ----------
        Paramaters:
            Seq(Bio.Seq)
            K(int)
        Returns:
            set(str)
    '''
    
    kmers = set()
    for i in range(len(Seq) - k +1):
        kmers.add(str(Seq[i:i+k]))
    return kmers

def count_occurrences(Seq,w):
    
    '''
        Counts the number of occurences of w in Seq
        --------
        Parameters:
            Seq (Bio.Seq)
            w (str)
        Returns:
            int : number of w in Seq
    '''
    
    count = 0
    for i in range(len(Seq) -len(w) +1):
        for j in range(len(w)):
            if Seq[i+j] != w[j]:
                break
        else:
            count += 1
    return count

def print_kmers_type (Seq, k):
    
    '''
        Prints the types of kmers in a sequence
        -----------
        Parameters:
        Seq(Bio.Seq)
        k(int)
    '''
    
    words=set()
    kmers = get_kmers(Seq,k)
    words = get_list_words('',k,set(),Seq)
    for kmer in words:
        if count_occurrences(Seq,kmer)==0:
            print(kmer + "is a nullomer")
    for kmer in kmers:
        m = count_occurrences(Seq,kmer)
        if m == 1:
            print(kmer + " is an hapax")
        else:
            if m == 2:
                print(kmer + " is a duplex")
            else:
                print(kmer + " is a repeat")

def print_nullomers (Seq,k):
    
    '''
        Prints the nullomers in a sequence
        -----------
        Parameters:
        Seq(Bio.Seq)
        k(int)
    '''
    
    words = get_list_words('',k,set())
    for kmer in words:
        if kmer not in Seq:
            print(kmer + "is a nullomer")

def print_hapaxes(Seq,k):
    
    '''
        Prints the hapaxes in a sequence
        -----------
        Parameters:
        Seq(Bio.Seq)
        k(int)
    '''
    
    kmers = get_kmers(Seq,k)
    for kmer in kmers:
        m = count_occurrences(Seq,kmer)
        if m == 1:
            print(kmer + " is an hapax")

def print_duplexes (Seq,k):
    
    '''
        Prints the duplexes in a sequence
        -----------
        Parameters:
        Seq(Bio.Seq)
        k(int)
    '''
    
    kmers = get_kmers(Seq,k)
    for kmer in kmers:
        m = count_occurrences(Seq,kmer)
        if m == 2:
            print(kmer + " is a duplex")
            
def print_repeats (Seq,k):
    
    '''
        Prints the repeats in a sequence
        -----------
        Parameters:
        Seq(Bio.Seq)
        k(int)
    '''
    
    kmers = get_kmers(Seq,k)
    for kmer in kmers:
        m = count_occurrences(Seq,kmer)
        if m > 2:
            print(kmer + " is a repeat")
                
def get_hapaxes(Seq,k):
    
    '''
      Returns a list of hapaxes in Seq
      -----------
      Parameters:
          Seq(Bio.Seq)
          k(int)
      --------
      Returns:
          sorted(list[str]) : hapaxes
      '''
    
    kmers = get_kmers(Seq,k)
    hapaxes = []
    for kmer in kmers:
        m = count_occurrences(Seq,kmer)
        if m == 1:
            hapaxes.append(kmer)
    return sorted(hapaxes)
def get_nullomers(Seq,k):
    
    '''
      Returns a list of nullomers in Seq
      -----------
      Parameters:
          Seq(Bio.Seq)
          k(int)
      --------
      Returns:
          sorted(list[str]) : nullomers
            
      '''
    
    nullomers = []
    words=get_list_words('',k,set())
    for kmer in words:
        if kmer not in Seq:
            nullomers.append(kmer)
    return sorted(nullomers)

def get_repeats(Seq,k):
    
    '''
      Returns a list of repeats in Seq
      -----------
      Parameters:
          Seq(Bio.Seq)
          k(int)
      --------
      Returns:
          sorted(list[str]) :repeats

    '''
    
    kmers = get_kmers(Seq,k)
    repeats = []
    for kmer in kmers:
        m = count_occurrences(Seq,kmer)
        if m>2:
            repeats.append(kmer)
    return sorted(repeats)
         
def get_duplexes(Seq,k):
    '''
       Returns a list of duplexes in a Seq
       ----------
       Parameters:
           Seq(Bio.Seq)
           k(int)
        --------
        Returns:
        sorted(list[str]) : duplexes
    '''
    kmers = get_kmers(Seq,k)
    duplexes = []
    for kmer in kmers:
        m = count_occurrences(Seq,kmer)
        if m==2:
            duplexes.append(kmer)
    return sorted(duplexes)
def get_kmers_type (Seq,k):
    '''
        Returns a list that contains three sublists (nullomers, hapaxes & repeats), also repeats contains a sublist(duplexes)
        -----------
        Parameters:
            Seq(Bio.Seq)
            k(int)
        --------
        Returns:
            list[list[str],list[str],list[str,list[str]]] : nullomers,hapaxes, repeats and duplexes
    '''
    
    kmers = get_kmers(Seq,k)
    words=set()
    kmer_type = []
    hapaxes = []
    repeats = []
    duplexes = []
    nullomers = []
    words=get_list_words('',k,set(),Seq)
    for kmer in words:
        if kmer not in Seq:
            nullomers.append(kmer)
    for kmer in kmers:
        m = count_occurrences(Seq,kmer)
        if m == 1:
            hapaxes.append(kmer)
        else:
            
            if m == 2:
                duplexes.append(kmer)
            else:    
                repeats.append(kmer)
    kmer_type.append(sorted(nullomers))
    kmer_type.append(sorted(hapaxes))
    repeats.append(sorted(duplexes))
    kmer_type.append(repeats)
    return kmer_type 

def print_list_words(prefix, k,Seq):
    '''
        Prints a set of words of length k in alphabet A,C,G,N,T
        -----------
        Parameters:
            prefix('')
            k(int)
            Seq(Bio.Seq)
    '''
    
    if len(prefix) == k:
        print(prefix)
    else:
        for a in get_alphabet(Seq):
            print_list_words(prefix + a, k,Seq)
def get_list_words(prefix, k, words_list,Seq):
    
    '''
        Returns a set of words of length k for an alphabet of A,C,G,N,T 
        -----------
        Parameters:
            prefix('')
            k(int)
            words_list(set())
            Seq(Bio.Seq)
        --------
        Returns:
            sorted(set(str)) : words of length k for an alphabet A,C,G,N,T
    '''
    
    alph=get_alphabet(Seq)
    if len(prefix)==0:
        words_list=set()
    if len(prefix) == k:
        words_list.add(prefix)
    else:
        for a in alph:
            get_list_words(prefix + a, k,words_list,Seq)
    return sorted(words_list)

def countN(Seq,char):
    
    '''
        Counts nucleotide in a Seq objects
        -----------
        Parameters:
            Seq(Bio.Seq)
            char("char")
        --------
        Returns:
            int : number of nucleotides in a Seq
    '''
    
    count = 0 
    for a in Seq:
        if a == char:
            count += 1
    return count

def get_positions(Seq,w):
    """
        Returns the starting postions in a reference sequence Seq where a word w occurs
        -----------
        Parameters:
            Seq(Bio.Seq) : the reference sequence
            w (str) : the searched word
        --------
        Returns:
            list[int] : the positions
    """
    positions = list()
    for i in range(len(Seq)):
        if Seq[i:i+len(w)] == w:
            positions.append(i)
    return positions

def dictionary_coverage_position(Seq,dictionary):
    
    '''
        Returns the positions coverage in a sequence Seq of a dictionary
        -----------
        Parameters:
            Seq(Bio.Seq) : the reference sequence
            dictionary[str] : dictionary of words (str)
        --------
        Returns:
            list[int] : dictionary coverage
    
    '''
    
    coverage = [0 for _ in range(len(Seq)) ]
    for w in dictionary:
        for pos in get_positions(Seq,w):
            for i in range(len(w)):
                coverage[pos + i] += 1
    return coverage
def coverage(Seq, dictionary):
    
    '''
        Returns the coverage of a Sequence
        ---------
        Parameters:
            Seq(Bio.Seq) : the reference sequence
            dictionary(str) :dictionary of words(str)
        --------
        Returns:
            int : coverage    
    '''
    
    coverage = dictionary_coverage_position(Seq,dictionary)
    coverage=((len(coverage) - coverage.count(0))  / len(coverage)) 
    return coverage

def mrl(Seq):
    
    """
        Calculates the maximal repeat length of a sequence Seq
        -----------
        Parameters:
            Seq(Bio.Seq) : the reference sequence
        -------
        Returns:
            int : maximal repeat length
            str : maximal repeat length kmer
            int : multiplicity
    """
    
    k = 0
    mrl = 0
    kmer_mrl = ''
    mult_mrl = 0
    next_k = True
    while next_k:
        #print(k, end='', sep=' ')
        k += 1
        next_k = False
        for kmer in get_kmers(Seq,k):
            mult = count_occurrences(Seq,kmer)
            if mult > 1:
                mrl = k
                kmer_mrl = kmer
                mult_mrl = mult
                next_k = True
    return mrl, kmer_mrl, mult_mrl

def mhl(Seq):
    
    """
        Calculate the minimal hapax length of a sequence Seq 
        ------------
        Parameters:
            Seq(Bio.Seq) : the reference sequence
        -------
        Returns:
        int: minimal hapax length
        str: minimal hapax length kmer
        int: multiplicity
        
    """
    
    k = 0
    mhl = 0
    kmer_mhl = ''
    mult_mhl = 0
    next_k = True
    while next_k:
        k += 1
        for kmer in get_kmers(Seq,k):
            mult = count_occurrences(Seq,kmer)
            if mult == 1:
                mhl = k
                kmer_mhl = kmer
                mult_mhl = mult
                next_k = False
    return mhl, kmer_mhl, mult_mhl

def mfl(Seq, alphabet = None):
    
    """
    Calculates the minimal forbidden length of a string s
    ----------
    Parameters:
        Seq(Bio.Seq): the reference sequence
        alphabet if passed
    
    """
    
    if alphabet == None:
        a = len(get_alphabet(Seq))
    else:
        a = len(alphabet)
        
    k = 0
    while True:
        k += 1
        kmers = get_kmers(Seq,k)
        if len(kmers) != a**k:
            return k

def get_multiplicity_distribution(Seq,k):
    """
        Returns the word multiplciity distribution of k-mers occurrig in the string s
         ------
        Parameters:
            s (str) : the input string
            k (int) : the length of the k-mers 
        -------
        Returns:
            dict[str,int] : a dictionary which associates multiplicity values to the k-mers in s
    """
    
    WMD = dict()
    for i in range(len(Seq) -k +1):
        w = str(Seq[i:i+k])
        WMD[w] = WMD.get(w,0) + 1     
    return WMD

def k_entropy(Seq,k):
    
    """
        Calculates the empirical entropy at word length k of a sequence Seq
        -----------
        Parameters:
            Seq(Bio.Seq) : the reference sequence
            k(int) : input k
        ----------
        Returns:
            int : input k
            int : - empirical entropy value 
    
    """
    
    distr = get_multiplicity_distribution(Seq,k)
    t = sum(distr.values())
    e = 0.0
    for v in distr.values():
        e += math.log(v / t, 2)
   
    return k,-e

def average_coverage (Seq,dictionary) :
    
    '''
        Returns the average positional coverage 
        ----------
        Parameters:
            Seq(Bio.Seq) : the reference sequence 
            dictionary(str) :dictionary of words(str)
        ----------
        Returns: 
            int : average positional coverage
    '''
    
    coverage=dictionary_coverage_position(Seq,dictionary)
    return mean(coverage)

def standard_deviation_coverage (Seq,dictionary) :
    
    '''
        Returns the standard deviation of positional coverage
        ----------
        Parameters:
            Seq(Bio.Seq) : the reference sequence 
            dictionary(str) :dictionary of words(str)
        ----------
        Returns int : standard deviation of positional coverage
    '''
    
    coverage=dictionary_coverage_position(Seq,dictionary)
    return stdev(coverage)


def covered_average_coverage (Seq,dictionary) :
    
    '''
        Returns the average positional coverage of covered positions
        ----------
        Parameters:
            Seq(Bio.Seq) : the reference sequence 
            dictionary(str) :dictionary of words(str)
        ----------
        Returns int : average positional coverage of covered positions
    '''
    
    coverage=dictionary_coverage_position(Seq,dictionary)
    return mean([ i for i in coverage if i > 0])

def covered_standard_deviation_coverage (Seq,dictionary) :
    
    '''
        Returns the standard deviation of positional coverage of covered positions
        ----------
        Parameters:
            Seq(Bio.Seq) : the reference sequence 
            dictionary(str) :dictionary of words(str)
        ----------
        Returns int : standard deviation of positional coverage of covered positions
    '''
    
    coverage=dictionary_coverage_position(Seq,dictionary)
    return stdev([ i for i in coverage if i > 0])

def wld(Seq, k_start, k_end):
    
    """
        Calculates the word lenght distribution of the sequence Seq for the given range of values of word length k
        ----
        Paramters:
            Seq (Bio.Seq) : the input string
s            k_start (int) : the initial word length
            k_end (int) : the final word length
        ----
        Returns:
            dict[int,int] : a dictionary which associates word lengths (key) to the number of k-mers at each length (value)
    """
    
    wld = dict()
    for k in range(k_start,k_end):
        wld[k] = len(get_kmers(Seq,k))
    return wld

def amd(Seq, k_start, k_end):
    
    """
        Calculates the average multiplcity distribution of the sequence Seq for the given range of values of word length k
        ----
        Paramters:
            Seq (Bio.Seq) : the input string
            k_start (int) : the initial word length
            k_end (int) : the final word length
        ----
        Returns:
            dict[int,int] : a dictionary which associates word lengths (key) to the average multiplicity at the specific word length (value)
    """
    
    if ((k_end-k_start)>len(Seq)):
        print("error k_end-k_start > Seq length");
        return null
    amd = dict()
    for k in range(k_start,k_end):
        amd[k] = statistics.mean( get_multiplicity_distribution(Seq,k).values() )
    return amd

def eed(Seq, k_start, k_end):
    """
        Calculates the empirical entropy distribution of the sequence Seq for the given range of values of word length k
        ----
        Paramters:
            Seq (Seq.Bio) : the input string
            k_start (int) : the initial word length
            k_end (int) : the final word length
        ----
        Returns:
            dict[int,int] : a dictionary which associates word lengths (key) to the empirical entropy at the specific word length (value)
    """
    
    eed = dict()
    for k in range(k_start,k_end):
        eed[k] = k_entropy(Seq,k)
    return eed

def wcmd(Seq,k):
    
    """
        Calculates the word co-multiplicity distribution of the sequence Seq for the given value of word length k
        ----
        Paramters:
            Seq (Bio.Seq) : the input string
            k (int) : the word length
        ----
        Returns:
            dict[int,int] : a dictionary which associates a multiplicity value (key) to the number of k-mers having such multiplicity (value)
    """
    
    distr = dict()
    mdistr = get_multiplicity_distribution(Seq,k)
    for m in mdistr.values():
        distr[m]= distr.get(m,0) + 1
    return distr

"""
   ______________________________________________________________________________________________
   W: THIS FUNCTION IS NOT OPTIMIZED. FOR SUFFIX ARRAYS WE USE PYSUFFIXARRAY WHICH IS WAY FASTER.

 def get_suffix_array(Seq):
    
    Construct the suffix array of the sequence Seq.
    
    pairs = list()
    for i in range(len(Seq)):
        pairs.append( (Seq[i:],i) ) 
    sa = list()
    for p in sorted(pairs):
        sa.append(p[1])
    return sa

"""

def longest_prefix_length(Seq, i, j):
    
    """
        Calculates the length of the longest common prefix between two suffixes, 
        the one in position i and the other in position j, of Seq.
        ----------------
        Parameters:
           Seq(Bio.Seq): the input string
           int i
           int j
        ---------------
        Returns:
           l: longest prefix between two suffixes, the one in position i and the other in position j
        
    """
    
    l = 0
    while (i+l < len(Seq)) and (j+l < len(Seq)):
        if Seq[i+l] != Seq[j+l]:
            break
        l += 1
    return l


def get_lcp(Seq):
    
    """
        Constructs the LCP array associated to the suffix array (sa) of the sequence Seq.
        The LCP value of the first suffix is set to be 0.
        ---------------------------
        Parameters:
            Seq(Bio.Seq): the input string
        -------------------------
        Returns:
           lcp: longest common prefix
    """
    
    sa=SuffixArray(Seq)
    lcp = list()
    lcp.append(0)
    for i in range(1,len(sa)):
        lcp.append( longest_prefix_length(Seq, sa[i], sa[i-1]) )
    return lcp

def print_sa_lcp(Seq):
    
    """
        Prints of the ordered suffixes together with the suffix array
        ----------------
        Parameters:
            Seq(Bio.Seq): the input string
    """
    
    sa=SuffixArray(Seq)
    lcp=get_lcp(Seq)
    print('index', 'suffixes' + ' '*(len(Seq)-len('suffixes')), 'SA', 'LCP', sep='\t')
    print('-'*45)
    for i in range(len(sa)):
        print(i, Seq[sa[i]:] + ' '*(sa[i]), sa[i], lcp[i], sep='\t')

def print_sa_lcp_region(Seq):
    
    """

        Prints the suffix array and the longest common prefix in a range i-j
        -----------
        Parameters: 
            Seq(Bio.Seq):the input string
    """
    
    sa=SuffixArray(Seq)
    lcp=get_lcp(Seq)

    i=int(input("Enter i: "))
    j=int(input("Enter j: "))
    
    print('-'*40)
    print('interval')
    for x in range(i,j):
        print(x,Seq[sa[x]:] +' '*(sa[x]), sa[x], lcp[x], sep='\t')
    #print('.'*40)


def distance_to_n(Seq,i):
    """

        Returns the distance of the first occurrency of the character 'N' from the index i in the string s
        -------
        Parameters:
            Seq(Bio.Seq):the input string
            i (int) : the start index
        -------
        Returns:
           j-i: the distance of the first occurrency of the...
    """
    
    j = i
    while (j<len(Seq)) and (Seq[j] != 'N'):
        j += 1
    return j - i

def get_ns_array(Seq):
    
    """
        Calls distance_to_n to count the distance in all the suffixes 
        -------
        Parameters:
            Seq(Bio.Seq):the input string
        -------
        Returns:
            the call to distance_to_n
    """
    
    sa = SuffixArray(Seq)
    return [ distance_to_n(Seq,sa[i]) for i in range(len(Seq)) ]


def print_sa_lcp_ns(Seq):
    
    """

        Prints the suffix array, the longest common prefix and nelsa
        -------
        Parameters:
           Seq(Bio.Seq):the input string
        -------
    """
    
    sa = SuffixArray(Seq)
    lcp = get_lcp(Seq)
    ns = get_ns_array(Seq)
        
    print('index', 'suffixes' + ' '*(len(Seq)-len('suffixes')), 'SA', 'LCP', 'NS', sep='\t')
    print('-'*60)
    for i in range(len(Seq)):
        print(i, Seq[sa[i]:] + ' '*(sa[i]), sa[i], lcp[i], ns[i], sep='\t')

def print_sa_lcp_ns_region(Seq,i,j):
    
    """

        Prints the suffix array, the longest common prefix and nelsa in a range i-j
        -------
        Parameters:
        Seq(Bio.Seq): the input string
            i (int) : the start index
            j (int) : the end index
    """
    
    sa = SuffixArray(Seq)
    lcp = get_lcp(Seq)
    ns = get_ns_array(Seq)
    
    print('-'*60)
    print('interval')
    for x in range(i,j):
        print(x,Seq[sa[x]:] +' '*(sa[x]), sa[x], lcp[x], ns[x], sep='\t')    

def fast_get_ns_array(Seq): 
    """

        Prints the suffix array, the longest common prefix and nelsa in a range i-j
    -------
    Parameters:
        Seq(Bio.Seq):the input string
    -------
    Returns: 
        ns: the list of nelsa
    """
    sa = SuffixArray(Seq)
    inv_sa = [0 for _ in range(len(sa))]
    for i in range(len(sa)):
        inv_sa[ sa[i] ] = i
    
    ns = [0 for _ in range(len(sa))]
    lastn = len(Seq)
    for i in range(len(Seq)-1,-1,-1):
        if Seq[i] == 'N':
            lastn = i
        ns[ inv_sa[i] ] = lastn - i
    return ns

class ESAIterator:
    __Seq = None
    __k = 0
    __sa = None
    __lcp = None
    __i = 0
    __j = 0
    
    def __init__(self, Seq, k, sa = None, lcp = None):
        self.__Seq = Seq
        self.__k = k
        
        if sa == None:
            self.build_sa()
        else:
            self.__sa = sa
            
        if lcp == None:
            self.build_lcp()
        else:
            self.__lcp = lcp

    def build_sa(self):
        print("building suffix array...")
        suffixes = list()
        for i in range(len(self.__Seq)):
            suffixes.append( (self.__Seq[i:] + self.__Seq[:i] ,i) ) 
        self.__sa = list()
        for suff in sorted(suffixes):
            self.__sa.append(suff[1])
        print('done')
        
        
    def longest_prefix_length(Seq, i, j):
        l = 0
        while (i+l < len(Seq)) and (j+l < len(Seq)):
            if Seq[i+l] != Seq[j+l]:
                break
            l += 1
        return l

    def build_lcp(self):
        print('building lcp array...')
        self.__lcp = list()
        self.__lcp.append(0)
        for i in range(1,len(self.__sa)):
            self.__lcp.append( ESAIterator.longest_prefix_length(self.__Seq, self.__sa[i], self.__sa[i-1]) )
        print('done')
    
    def get_sa(self):
        return self.__sa
    
    def get_lcp(self):
        return self.__lcp
        
    def __iter__(self):
        return self
    def __next__(self):
        if self.__i < len(self.__Seq):
            self.__i = self.__j
            
            while (self.__i < len(self.__Seq)) and  (self.__sa[self.__i] > len(self.__Seq) - self.__k - 1):
                self.__i += 1
            if self.__i == len(self.__Seq):
                raise StopIteration
            self.__j = self.__i+1
            while ( self.__j < len(self.__Seq) ) and (self.__lcp[self.__j] >= self.__k):
                self.__j += 1
            ret = self.__Seq[ self.__sa[self.__i] : self.__sa[self.__i] + self.__k ]
            #self.__i = self.__j #!!!!!!
            return ret
        else:
            raise StopIteration
            
    def multiplicity(self):
        return self.__j - self.__i
    
    def positions(self):
        return self.__sa[self.__i : self.__j]

class NESAIterator:
    __Seq = None
    __k = 0
    __sa = None
    __lcp = None
    __ns = None
    __i = 0
    __j = 0
    
    def __init__(self, Seq, k, sa = None, lcp = None, ns = None):
        self.__Seq = Seq
        self.__k = k
        
        if sa == None:
            self.build_sa()
        else:
            self.__sa = sa
            
        if lcp == None:
            self.build_lcp()
        else:
            self.__lcp = lcp
            
        if ns == None:
            self.build_ns()
        else:
            self.__ns = ns

    def get_k(self):
        return self.__k
    
    def reset(self):  
        self.__i = 0
        self.__j = 0
            
    def build_sa(self):
        print("building suffix array...")
        suffixes = list()
        for i in range(len(self.__Seq)):
            suffixes.append( (self.__Seq[i:] + self.__Seq[:i] ,i) ) 
        self.__sa = list()
        for suff in sorted(suffixes):
            self.__sa.append(suff[1])
        print('done')
        
        
    def longest_prefix_length(Seq, i, j):
        l = 0
        while (i+l < len(Seq)) and (j+l < len(Seq)):
            if Seq[i+l] != Seq[j+l]:
                break
            l += 1
        return l

    def build_lcp(self):
        print('building lcp array...')
        self.__lcp = list()
        self.__lcp.append(0)
        for i in range(1,len(self.__sa)):
            self.__lcp.append( NESAIterator.longest_prefix_length(self.__Seq, self.__sa[i], self.__sa[i-1]) )
        print('done')
    
    def build_ns(self):
        print('building ns array...')
        inv_sa = [0 for _ in range(len( self.__sa))]
        for i in range(len(self.__sa)):
            inv_sa[  self.__sa[i] ] = i

        self.__ns = [0 for _ in range(len( self.__sa))]
        lastn = len(self.__Seq)
        for i in range(len(self.__Seq)-1,-1,-1):
            if self.__Seq[i] == 'N':
                lastn = i
            self.__ns[ inv_sa[i] ] = lastn - i
        print('done')
        
    def get_sa(self):
        return self.__sa
    
    def get_lcp(self):
        return self.__lcp
    
    def get_ns(self):
        return self.__ns
        
        
    def __iter__(self):
        return self
    def __next__(self):
        if self.__i < len(self.__Seq):
            self.__i = self.__j
            
            while (self.__i < len(self.__Seq)) and  ( (self.__sa[self.__i] > len(self.__Seq) - self.__k - 1) or (self.__ns[self.__i] < self.__k) ):
                self.__i += 1
            if self.__i == len(self.__Seq):
                raise StopIteration
            self.__j = self.__i+1
            while ( self.__j < len(self.__Seq) ) and (self.__lcp[self.__j] >= self.__k) and (self.__ns[self.__i] >= self.__k) :
                self.__j += 1
            ret = self.__Seq[ self.__sa[self.__i] : self.__sa[self.__i] + self.__k ]
            #self.__i = self.__j #!!!!!!
            return ret
        else:
            raise StopIteration
            
    def multiplicity(self):
        return self.__j - self.__i
    
    def positions(self):
        return self.__sa[self.__i : self.__j]
    
def RDD(Seq,w):
    
    """
        Extracts the recurrence distance ditribution (RDD) of the word w in Sequence.
        Given the starting postions of two occurences of w, p1 and p2, the reucrrence distance is calculated as
        p1 - p2
        such that consecutive occurrences are at distance 1.
        ----
        Parameters:
            Seq(Bio.Seq):the input string
            w (str) : the searched substring
        ----
        Returns:
            dict[int,int] : a dictionary mapping recurrence distances to the number of times they occur
    """
    
    if(w in Seq):
        pos = sorted(get_positions(Seq,w))
        rdd = dict()
        for i in range(2,len(pos)):
            l = pos[i] - pos[i-1] 
            rdd[l] = rdd.get(l,0) + 1
        return rdd
    else:
        print('w not contained in sequence Seq')

def plot_RDD(rdd,title):
    
    """
        Plot an RDD distribution adding missing recurring distances, 
        between 1 and the original maximum distance,
        by adding a value of zero in correspondence of them.
        -------------
        Parameters:
            rdd: recurrence distance distribution
            title: title of the plot
    """
            
    # if a value equal to zero for the missing recurrence distances
    for d in range(0,max(rdd.keys())):
        rdd[d] = rdd.get(d,0) + 0
        
    # module can be imported by uring aliases to refer them
    import matplotlib.pyplot as plt 
    
    # set the figure size
    plt.rcParams['figure.figsize'] = [20, 6]
    # assign height of bars
    bar_values = [v for k,v in sorted(rdd.items())]
    # plot with specific values on the x-axis that are associated to the height
    plt.bar(sorted(rdd.keys()), bar_values, width=1.0)
    # set the label on the y-axis
    plt.ylabel('Number of pairs') 
    # set the label on the x-axis
    plt.xlabel('Recurrence distance')
    # set a title for the chart
    plt.title(title)
    # plot the chart
    plt.show()

    
def aRDD(Seq,k):
    
    """
        Computer the average recurrence distance distribution of the complete set of k-mers occuring in s.
        ----
        Parameters:
            Seq(Bio.Seq):the input sequence 
            k (int) : the word length of the k-mers for which extract the RDD
        ----
        Returns:
            dict[int,float] : a dictionary mapping recurrence distances to the average number of times they occur
    """
    
    ardd = dict()
    kmers = get_kmers(Seq,k)
    for kmer in kmers:
        rdd = RDD(Seq,kmer)
        for distance,value in rdd.items():
            ardd[distance] = ardd.get(distance,0) + value
    for d,v in ardd.items():
        ardd[d] = ardd[d] / len(kmers)
    return ardd

def genome_file_info (ifileFASTA, ifileGFF3):
    
    """
       Returns the sequence coverage,total length,non-coding length and protein-coding length of the sequence.
       ------------
       Parameters:
           ifileFASTA: the .fna file
           ifileGFF3: the .gff3 file
    """
    
    genome = ''
    for line in open(ifileFASTA, 'r'):
        if line.strip()[0] != '>':
            genome += line.strip()

    coverage = [0 for i in range(len(genome))]
        
    for line in open(ifileGFF3, 'r'):
        if line[0] != "#":
            cc = line.split('\t')
            if len(cc) >= 6:
                if (cc[2] == 'gene'):# and (cc[6] == '+'): # we calculate the coverage of both strands as a single strand
                    start = int(cc[3])
                    end = int(cc[4])
                    for i in range(start, end):
                        coverage[i] += 1

    print('sequence coverage',   (len(coverage) - coverage.count(0)) / len(coverage))

# sequence of non-coding portion of the genome
    ncseq = ''.join([ genome[i] for i in range(len(genome)) if coverage[i] == 0 ])

# sequence of coding portion of the genome
    cseq = ''.join([ genome[i] for i in range(len(genome)) if coverage[i] > 0 ])

    print('total length', len(genome),', non-coding length',len(ncseq), ', protein-coding length', len(cseq))


def list_words_2(prefix, k, words):

    '''
        Returns a set of words of length k for an alphabet of A,C,G,N,T 
        -----------
        Parameters:
            prefix('')
            k(int)
            words(set())
        -----------
    '''
    
    if len(prefix) == k:
        words.add(prefix)
    else:
        for a in nuc_alphabet():
            list_words_2(prefix + a, k, words)
            
def plot_wmd (k,Seq):

    """
       Plots word multiplicity distance.
       ----------
       Parameters
           k(int)
           Seq(Bio.Seq):the input sequence
    """
    
    kmers = set()
    list_words_2('',k, kmers)
    kmers = sorted(kmers)

    ncseq_wmd = dict()
    for kmer in kmers:
        ncseq_wmd[kmer] = count_occurrences(Seq, kmer)


    bar_values = [v for k,v in sorted(ncseq_wmd.items())]
    plt.rcParams['figure.figsize'] = [20, 6]
    plt.bar(kmers, bar_values)
    plt.xticks(kmers, kmers)
    plt.ylabel('multiplicity')
    plt.xlabel('words')
    plt.title('Word multiplicity distribution')
    plt.show()

    
def open_genome_fasta(filefasta):
    
    """
       Returns genome sequence in format Bio.Seq
       -----------
       Parameters:
           filefasta: The fasta file
    """
    
    for record in SeqIO.parse(filefasta, "fasta"):
        my_genome = record.seq
    return my_genome
