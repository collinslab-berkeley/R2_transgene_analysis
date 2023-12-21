import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from itertools import product
import pysam, os, sys, pickle
from dataclasses import dataclass
from copy import deepcopy
from Bio.Seq import Seq
from Bio import SeqIO
from matplotlib.gridspec import GridSpec
from scipy.signal import savgol_filter
from fuzzysearch import find_near_matches
from matplotlib import patches
from scipy import stats
from process_reads_v3 import *

@dataclass
class read_part:
    cigar: str
    read_start: int
    read_end: int
    chrom: str=''
    strand: chr=''
    ref_start: int = -1
    ref_end: int = -1
    mapq: int = -1
    
    def get_length(self):
        return self.read_end-self.read_start

pd.set_option('mode.chained_assignment', None)

def is_rDNA(chrom, coord, rDNA_ref='KY962518.1'):
    is_rDNA = chrom == rDNA_ref
    chr13 = chrom == 'chr13' and 5770548 < coord < 9348041
    chr14 = chrom == 'chr14' and 2099537 < coord < 2817811
    chr15 = chrom == 'chr15' and 2506442 < coord < 4707485
    chr21 = chrom == 'chr21' and 3108298 < coord < 5612715
    chr22 = chrom == 'chr22' and 4793794 < coord < 5720650
    is_rDNA = is_rDNA or chr13 or chr14 or chr15 or chr21 or chr22
    return is_rDNA

def get_bestmatch(matches):
    if len(matches) == 0:
        return 0
    elif len(matches) == 1:
        return matches[0]
    else:
        LD_list = np.array([match.dist for match in matches])
        # is there single best match?
        if np.sum(LD_list == LD_list.min()) == 1: 
            return matches[np.argmin(LD_list)]
        else:
            return 0
    
def rc(seq):
    return str(Seq(seq).reverse_complement())

def is_polyA(seq):
    base_dict = Counter(seq)
    return base_dict['A']/len(seq) > 0.8 and len(seq)>10

def is_polyG(seq):
    base_dict = Counter(seq)
    return base_dict['G']/len(seq) > 0.8 and len(seq)>10

def get_transgene_coverage(dict_in, transgene_refname, transgene_len, flank_len=840):
    pos_list = np.array([], dtype=int)
    end = transgene_len + flank_len

    for read_ID in dict_in:
        temp_list = np.array([], dtype=int)
        for read in range(2):
            parts = [x for x in dict_in[read_ID][read].keys() if 'pt' in x]
            for part in parts:
                read_part = dict_in[read_ID][read][part]
                if read_part.chrom == transgene_refname:
                    temp_list = np.append(temp_list,np.arange(read_part.ref_start,read_part.ref_end))

        pos_list = np.append(pos_list,np.unique(temp_list))
        
    bins = np.arange(flank_len, end+1)
    y = plt.hist(pos_list, bins=bins)
    plt.close()
    
    return y[0][:-1]

def classify_junctions(dict_in, transgene_refname, transgene_len, buffer=4, flank_len=840):
    end = transgene_len+flank_len
    
    # dictionaries. read_ID: (read_num, first_pt, second_pt)
    typeI = {}
    typeII = {}
    typeIIb = {}
    typeIIIa = {} # match-unmapped
    typeIIIb = {} # unmapped-match

    # dictionaries. read_ID: (read_num, part) 
    typeIVa = {} # 5' junctions
    typeIVb = {} # 3' junctions
    
    for read_ID in dict_in:
        temp_list = np.array([], dtype=int)
        for read in range(2):
            cigar_tuple = [dict_in[read_ID][read][x].cigar for x in dict_in[read_ID][read].keys() if 'pt' in x]
            # find match-match, match-unmapped, and unmapped-match junctions
            for i in range(len(cigar_tuple)-1):
                if dict_in[read_ID][read]['pt%i' % (i+1)].chrom == transgene_refname or dict_in[read_ID][read]['pt%i' % (i+2)].chrom == transgene_refname:
                    if cigar_tuple[i] == 'match' and cigar_tuple[i+1] == 'match':
                        typeI[read_ID] = (read,i,i+1)
                    if cigar_tuple[i] == 'match' and cigar_tuple[i+1] == 'unmapped' and i+3>len(cigar_tuple):
                        typeIIIa[read_ID] = (read,i,i+1)
                    if i==0 and cigar_tuple[i] == 'unmapped' and cigar_tuple[i+1] == 'match':
                        typeIIIb[read_ID] = (read,i,i+1)
            # find match-unmapped-match junctions
            for j in range(len(cigar_tuple)-2):
                if dict_in[read_ID][read]['pt%i' % (j+1)].chrom == transgene_refname or dict_in[read_ID][read]['pt%i' % (j+3)].chrom == transgene_refname:
                    if cigar_tuple[j] == 'match' and cigar_tuple[j+1] == 'unmapped' and cigar_tuple[j+2] == 'match':
                        typeII[read_ID] = (read,j,j+2)
                    elif cigar_tuple[j] == 'match' and cigar_tuple[j+1] == 'insertion' and cigar_tuple[j+2] == 'match':
                        typeIIb[read_ID] = (read,j,j+2)

            # find junctions crossing transgene
            parts = [x for x in dict_in[read_ID][read].keys() if 'pt' in x]
            for part in parts:
                entry = dict_in[read_ID][read][part]
                if entry.chrom == transgene_refname:
                    if flank_len-buffer > entry.ref_start and flank_len+buffer < entry.ref_end:
                        typeIVa[read_ID] = (read,part)
                    elif end-buffer > entry.ref_start and end+buffer < entry.ref_end:
                        typeIVb[read_ID] = (read,part)
                        
    return typeI, typeII, typeIIb, typeIIIa, typeIIIb, typeIVa, typeIVb

def draw_connection(ax, start, end, flank_len=840):
    start -= flank_len
    end -= flank_len
    w = end-start
    xy = (start+w/2,0)
    ax.add_patch(patches.Arc(xy, width=w, height=w/1000))
    
def plot_gapped_insertions(ax, dict_in, junction_dict, reference, tg_length, flank_len=840):

    end = tg_length+flank_len
    coords = []

    for matepair in junction_dict:
        read, up, down = junction_dict[matepair]
        upstream = dict_in[matepair][read]['pt%i' % (up+1)]
        downstream = dict_in[matepair][read]['pt%i' % (down+1)]

        if upstream.chrom == reference and downstream.chrom == reference and upstream.ref_end > (flank_len+5) and downstream.ref_start < (end-10) and \
            (upstream.ref_end+14) < downstream.ref_start and upstream.strand == downstream.strand:
            coords.append((upstream.ref_end, downstream.ref_start))
            # print(matepair)
            
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_yticks([])
    ax.set_xlim(0,tg_length)
    ax.set_ylim(0,1.2)
    
    for coord in coords:
        draw_connection(ax,coord[0],coord[1])
        # print(coord[0],coord[1])
    ax.text(50,1,'$n=%i$' % len(coords), fontsize=12)
    
    return coords


def classify_typeI_typeII(junction_dict, dict_in, transgene_refname, transgene_len,
    flank_len=840, rDNA_ref='KY962518.1', target_loc=11679):
    end = transgene_len+flank_len
    same_strand = {}
    ss_other = []
    opposite_strand = {}
    os_other = []

    for matepair in junction_dict:
        read, up, down = junction_dict[matepair]
        upstream = dict_in[matepair][read]['pt%i' % (up+1)]
        downstream = dict_in[matepair][read]['pt%i' % (down+1)]

        if upstream.strand == downstream.strand:

            # gaps and indels
            if upstream.chrom == transgene_refname and downstream.chrom == transgene_refname and upstream.ref_end > (flank_len+5) and downstream.ref_start < (end-2) and \
                upstream.ref_end <= downstream.ref_start:
                if (upstream.ref_end+14) < downstream.ref_start:
                    same_strand[matepair] = ('gap',read,up,down)
                else:
                    same_strand[matepair] = ('indel',read,up,down)

            # 5' junction, rDNA upstream
            elif ((upstream.chrom == transgene_refname and (upstream.ref_end < flank_len or upstream.ref_start > end)) or is_rDNA(upstream.chrom,upstream.ref_end)) and \
                (downstream.chrom == transgene_refname and downstream.ref_start < (end-10) and downstream.ref_end > flank_len+3):
                same_strand[matepair] = ('5p_junction_rDNA',read,up,down)

            # 5' junction, other genomic upstream
            elif (downstream.chrom == transgene_refname and downstream.ref_start < (end-10) and downstream.ref_end > flank_len+5) and upstream.chrom != transgene_refname:
                seq = dict_in[matepair][read]['seq'][upstream.read_start:upstream.read_end]
                if not is_polyG(seq) and not is_polyA(seq):
                    same_strand[matepair] = ('5p_junction_genome',read,up,down)

            # internal joins and tandem insertions
            elif downstream.chrom == upstream.chrom:
                if upstream.ref_end > end-5:
                    same_strand[matepair] = ('tandem',read,up,down)
                else:
                    same_strand[matepair] = ('internal_join',read,up,down)

            # 3' junctions (rDNA)          
            elif (upstream.chrom == transgene_refname and upstream.ref_end < end+5) and is_rDNA(downstream.chrom,downstream.ref_start):
                same_strand[matepair] = ('3p_junction_rDNA',read,up,down)

            # 3' junctions (genome)          
            elif (upstream.chrom == transgene_refname and upstream.ref_end < end+5) and not is_rDNA(downstream.chrom,downstream.ref_start):
                same_strand[matepair] = ('3p_junction_genome',read,up,down)
                
            # other RNA inserted         
            elif upstream.chrom != transgene_refname and (is_rDNA(downstream.chrom,downstream.ref_start) or
                                                      (downstream.chrom==transgene_refname and downstream.ref_start>end-5)):
                seq = dict_in[matepair][read]['seq'][upstream.read_start:upstream.read_end]
                if not is_polyA(seq):
                    # same_strand[matepair] = ('targetsite_otherRNA',read,up,down)
                    print(matepair, 'other RNA inserted or repair artifact')
                                                      
            else:
                print(matepair,read)
                ss_other.append(matepair)

        else: # opposite strand

            # snap-backs
            if upstream.chrom == transgene_refname and downstream.chrom == transgene_refname:
                opposite_strand[matepair] = ('snapback_cDNA',read,up,down)

            # 3' junction, rDNA downstream (other than target site)
            elif upstream.chrom == transgene_refname and is_rDNA(downstream.chrom,downstream.ref_end):
                opposite_strand[matepair] = ('3p_junction_rDNA',read,up,down)

            # 3' junction, genomic sequence downstream (other than target site)
            elif upstream.chrom == transgene_refname and not is_rDNA(downstream.chrom,downstream.ref_end):
                opposite_strand[matepair] = ('3p_junction_genome',read,up,down)

            # rDNA snap-back
            elif downstream.chrom == transgene_refname and upstream.chrom == rDNA_ref and np.abs(upstream.ref_start-target_loc)<50:
                opposite_strand[matepair] = ('snapback_rDNA',read,up,down)

            # join to somewhere else in rDNA
            elif downstream.chrom == transgene_refname and is_rDNA(upstream.chrom,upstream.ref_end):
                opposite_strand[matepair] = ('5p_junction_rDNA',read,up,down)

            # join to somewhere else in genome
            elif downstream.chrom == transgene_refname and not is_rDNA(upstream.chrom,upstream.ref_end):
                opposite_strand[matepair] = ('5p_junction_genome',read,up,down)
                
            else:
                print(matepair,'other')
                os_other.append(matepair)
        
    return same_strand, ss_other, opposite_strand, os_other

def classify_5p_junctions(dict_in, junction_dict, flank_len=840, rDNA_ref='KY962518.1', target_loc=11679):
    _5p_labels = ['5p_junction_rDNA', '5p_junction_genome', 'snapback_cDNA', 'snapback_rDNA', 'tandem']
    classification_dict = {}
    rDNA_join_sites = {'FL':[], 'Trunc':[]}
    tg_start_sites = []    
    
    for matepair in junction_dict:
        label, read, up, down = junction_dict[matepair]
        upstream = dict_in[matepair][read]['pt%i' % (up+1)]
        downstream = dict_in[matepair][read]['pt%i' % (down+1)]
        if label in _5p_labels:
            tg_start_sites.append(np.max([downstream.ref_start, flank_len]))

            prefix = 'Trunc' # default is truncated
            if downstream.ref_start <= flank_len: # full-length
                prefix = 'FL'            

            if label == '5p_junction_rDNA':
                if prefix == 'FL':
                    rDNA_join_sites['FL'].append(upstream.ref_end)
                else:
                    rDNA_join_sites['Trunc'].append(upstream.ref_end)
                    
                if upstream.chrom == rDNA_ref: 
                    if target_loc - 4 <= upstream.ref_end <= target_loc + 2:
                        classification_dict[matepair] = '%s_join' % prefix
                    elif target_loc - 4 > upstream.ref_end:
                        classification_dict[matepair] = '%s_del' % prefix
                    else:
                        classification_dict[matepair] = '%s_dup' % prefix
                else:
                    print(matepair, upstream.chrom)
            elif label == '5p_junction_genome':
                classification_dict[matepair] = '%s_other' % prefix
            elif label == 'snapback_cDNA':
                classification_dict[matepair] = '%s_cDNA_snapback' % prefix
            elif label == 'snapback_rDNA':
                classification_dict[matepair] = '%s_rDNA_snapback' % prefix
            elif label == 'tandem':
                classification_dict[matepair] = '%s_tandem' % prefix
            else: # in case we missed a label <-- for debugging
                print(label)
    
    return classification_dict, rDNA_join_sites, tg_start_sites

def classify_typeIVb_junctions(dict_in, junction_dict, transgene_end, upstream_seq='TGTTCGG', downstream_seq='TAGCCAA'):
    IVS_dict = {}
    classification_dict = {}
    initiation_sites = []

    for read_ID in junction_dict:
        initiation_sites.append(transgene_end)
        classification_dict[read_ID] = 'on-target'
        read, part = junction_dict[read_ID]
        entry = dict_in[read_ID][read][part]
        seq = dict_in[read_ID][read]['seq'][entry.read_start:entry.read_end]
        
        # identify known upstream and downstream sequence to extract intervening sequence
        upstream_match = get_bestmatch(find_near_matches(upstream_seq, seq, max_l_dist=1))
        downstream_match = get_bestmatch(find_near_matches(downstream_seq, seq, max_l_dist=1))
        if upstream_match and downstream_match and upstream_match.end < downstream_match.start:
            if seq[upstream_match.end-1:downstream_match.start+4] == '':
                print(read_ID)
            IVS_dict[read_ID] = seq[upstream_match.end-1:downstream_match.start+4]
        elif downstream_match and downstream_match.start<11:
            upstream_match = find_near_matches('GAAA', seq[:downstream_match.start], max_l_dist=0)
            if upstream_match:
                if seq[upstream_match[0].start:downstream_match.start+4] == '':
                    print(read_ID)
                    pass
                IVS_dict[read_ID] = seq[upstream_match[0].start:downstream_match.start+4]
        elif upstream_match and len(seq)-upstream_match.end<11:
            downstream_match = find_near_matches('TAGCC', seq[upstream_match.end:], max_l_dist=0)
            if downstream_match:
                if seq[upstream_match.end-1:upstream_match.end-1+downstream_match[0].end] == '':
                    print(read_ID)
                    pass
                IVS_dict[read_ID] = seq[upstream_match.end-1:upstream_match.end-1+downstream_match[0].end]
        else:
            upstream_match = get_bestmatch(find_near_matches(rc(upstream_seq), seq, max_l_dist=1))
            downstream_match = get_bestmatch(find_near_matches(rc(downstream_seq), seq, max_l_dist=1))
            if upstream_match and downstream_match:
                IVS_dict[read_ID] = rc(seq[downstream_match.end-4:upstream_match.start+1])

    return classification_dict, initiation_sites, IVS_dict


def classify_3p_junctions(dict_in, junction_dict, transgene_len,
    flank_len=840, rDNA_ref='KY962518.1', target_loc=11679):
    end = transgene_len + flank_len
    upstream_seq = 'TGTTCGG'
    downstream_seq = 'TAGCCAA'
    _3p_labels = ['3p_junction_rDNA', '3p_junction_genome']
    classification_dict = {}
    initiation_temp = []
    IVS = {}
    uncertain = {}

    for matepair in junction_dict:
        label, read, up, down = junction_dict[matepair]
        upstream = dict_in[matepair][read]['pt%i' % (up+1)]
        downstream = dict_in[matepair][read]['pt%i' % (down+1)]
        if label in _3p_labels:
            seq = dict_in[matepair][read]['seq']
            internal_initiation = upstream.ref_end < end - 5 # false if started at 3' end of transcript
            if label == '3p_junction_rDNA':
                if internal_initiation:
                    if downstream.chrom==rDNA_ref and target_loc-5<downstream.ref_end<target_loc+5:
                        initiation_temp.append(upstream.ref_end)
                        classification_dict[matepair] = 'on-target internal'
                    else:
                        classification_dict[matepair] = 'uncertain'
                        uncertain[matepair] = junction_dict[matepair]
                else:
                    seq = dict_in[matepair][read]['seq']
                    # identify known upstream and downstream sequence to extract intervening sequence
                    upstream_match = get_bestmatch(find_near_matches(upstream_seq, seq, max_l_dist=1))
                    downstream_match = get_bestmatch(find_near_matches(downstream_seq, seq, max_l_dist=1))
                    if upstream_match and downstream_match:
                        IVS[matepair] = seq[upstream_match.end-1:downstream_match.start+4]
                        classification_dict[matepair] = 'on-target'
                    elif downstream_match and not upstream_match and downstream_match.start<11:
                        upstream_match = find_near_matches('GAAA', seq[:downstream_match.start], max_l_dist=0)
                        if upstream_match:
                            IVS[matepair] = seq[upstream_match[0].start:downstream_match.start+4]
                            classification_dict[matepair] = 'on-target'
                    elif upstream_match and not downstream_match and len(seq)-upstream_match.end<11:
                        downstream_match = find_near_matches('TAGCC', seq[upstream_match.end:], max_l_dist=0)
                        if downstream_match:
                            IVS[matepair] = seq[upstream_match.end-1:upstream_match.end-1+downstream_match[0].end]
                            classification_dict[matepair] = 'on-target'
                    else:
                        upstream_match = get_bestmatch(find_near_matches(rc(upstream_seq), seq, max_l_dist=1))
                        downstream_match = get_bestmatch(find_near_matches(rc(downstream_seq), seq, max_l_dist=1))
                        if upstream_match and downstream_match:
                            IVS[matepair] = rc(seq[downstream_match.end-4:upstream_match.start+1])
                            classification_dict[matepair] = 'on-target'
                        else:
                            classification_dict[matepair] = 'rDNA off-target'
            elif not internal_initiation:
                if is_polyA(seq[downstream.read_start:downstream.read_end]):
                    classification_dict[matepair] = 'indeterminate'
                else:
                    classification_dict[matepair] = 'genomic off-target'
            else:
                classification_dict[matepair] = 'uncertain'
                uncertain[matepair] = junction_dict[matepair]
                
    return classification_dict, initiation_temp, IVS, uncertain

def get_microhomology(dict_in, junction_dict, transgene_seq, rDNA_seq, transgene_refname, rDNA_refname='KY962518.1'):
    microhomology = {}
    upstream_nt = {}
    downstream_nt = {}
    join_chr = {}
    join_loc = {}
    ref_dict = {transgene_refname: transgene_seq, rDNA_refname: rDNA_seq}
    
    for matepair in junction_dict:
        read, up, down = junction_dict[matepair]
        upstream = dict_in[matepair][read]['pt%i' % (up+1)]
        downstream = dict_in[matepair][read]['pt%i' % (down+1)]
        
        if upstream.chrom in ref_dict and downstream.chrom in ref_dict:
            upstream_ref = ref_dict[upstream.chrom]
            downstream_ref = ref_dict[downstream.chrom]
            
            # check how many bp match downstream
            j = 0
            upstream_seq = upstream_ref[upstream.ref_end:]
            downstream_seq = downstream_ref[downstream.ref_start:]
            while upstream_seq[j]==downstream_seq[j]:
                j+=1
            downstream_nt[matepair] = downstream_seq[j:j+5]

            # check how many bp match upstream
            k = 0
            upstream_seq = upstream_ref[:upstream.ref_end]
            downstream_seq = downstream_ref[:downstream.ref_start]
            while downstream_seq[-1-k]==upstream_seq[-1-k]:
                k+=1

            microhomology[matepair] = j+k
            upstream_nt[matepair] = upstream_seq[-5-k:-1-k]+upstream_seq[-1-k]
            
            join_chr[matepair] = upstream.chrom
            join_loc[matepair] = upstream.ref_end-k-1

    return microhomology, upstream_nt, downstream_nt, join_chr, join_loc
