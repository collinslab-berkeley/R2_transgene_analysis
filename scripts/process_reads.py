import numpy as np
import pysam, os, pickle
from Bio.Seq import Seq
from Bio import SeqIO
from fuzzysearch import find_near_matches
import pysam, os, pickle
from dataclasses import dataclass
from copy import deepcopy

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

def get_orientation(reverse_bool):
    if reverse_bool:
        return '-'
    else:
        return '+'
    
def rc(seq):
    return str(Seq(seq).reverse_complement())

def get_alignments(seq_list, ref_file, fasta_name, sam_name):
    with open(fasta_name, 'w') as fa_out:
        for seq in seq_list:
            fa_out.write('>%s_%s\n' % (seq[0],seq[1]))
            fa_out.write(seq[2]+'\n')
            
    os.system('bwa mem -B 6 %s %s > %s' % (ref_file, fasta_name, sam_name))

    return pysam.AlignmentFile(sam_name, 'r')

def populate_read_dict(samfile_obj):

    all_reads = {}
    for read in samfile_obj.fetch():
        if read.query_name and not read.is_supplementary:
            if read.query_name in all_reads:
                all_reads[read.query_name].append(read)
            else:
                all_reads[read.query_name] = [read]

    read_dict = {}

    for matepair in all_reads:
        read_list = [None]*2
        for read in all_reads[matepair]:
            read_num = read.is_read2
            read_list[read_num] = {}
            read_list[read_num]['seq'] = read.seq
            idx = 1
            nt_count = 0
            nt_mapped = 0
            if read.cigar:
                for operation, length in read.cigar:
                    if operation in [0,1,4]:
                        read_list[read_num]['pt%i' % idx] = {} # create a dictionary for new part of read
                        if operation == 0: # matches
                            read_list[read_num]['pt%i' % idx] = read_part(cigar='match',
                                                                          read_start=nt_count,
                                                                          read_end=nt_count+length,
                                                                          chrom=read.reference_name,
                                                                          strand=get_orientation(read.is_reverse),
                                                                          ref_start=read.pos+nt_mapped,
                                                                          ref_end=read.pos+nt_mapped+length,
                                                                          mapq=read.mapq)
                            nt_count += length
                            nt_mapped += length
                        elif operation == 1: # insertion
                            read_list[read_num]['pt%i' % idx] = read_part(cigar='insertion',
                                                                          read_start=nt_count,
                                                                          read_end=nt_count+length)
                            nt_count += length
                        elif operation == 4: # soft-clipped
                            read_list[read_num]['pt%i' % idx] = read_part(cigar='unmapped',
                                                                          read_start=nt_count,
                                                                          read_end=nt_count+length)
                            nt_count += length
                        idx += 1 # move to next index in dict if portion of read is clipped, aligned, or an insertion
                    elif operation == 2: # deletion
                        nt_mapped += length
            else:
                read_list[read_num]['pt%i' % idx] = read_part(cigar='unmapped',
                                                              read_start=0,
                                                              read_end=len(read.seq))


        read_dict[matepair] = read_list
    
    return read_dict


def realign_reads(dict_in, ref_file, fasta_name='temp.fa', sam_name='temp.sam'):
    dict_copy = deepcopy(dict_in)

    seqs2write = []
    for matepair in dict_copy:
        for read in range(2):
            # if not fully aligned in one part
            if not (dict_copy[matepair][read]['pt1'].cigar == 'match' and 'pt2' not in dict_copy[matepair][read]):
                seqs2write.append([matepair,read,dict_copy[matepair][read]['seq']])

    realigned_reads = get_alignments(seqs2write, ref_file, fasta_name, sam_name)

    for read in realigned_reads.fetch():
        read_ID = read.query_name[:-2]
        read_num = int(read.query_name[-1])
        temp_dict = {}
        
        if read.cigar and not read.is_supplementary:
            temp_dict['seq'] = read.seq

            # if new portion is fully mapped, replace all info with new info
            if len(read.cigar) == 1:
                operation, length = read.cigartuples[0]
                temp_dict['pt1'] = read_part(cigar='match',
                                             read_start=0,
                                             read_end=length,
                                             chrom=read.reference_name,
                                             strand=get_orientation(read.is_reverse),
                                             ref_start=read.pos,
                                             ref_end=read.pos+length,
                                             mapq=read.mapq)


                # replace alignment in dict
                dict_copy[read_ID][read_num] = temp_dict

            # if new portion has 2 or more parts, populate 
            elif len(read.cigar) > 1:
                new_cigar = read.cigar
                new_bp_aligned = np.sum([x[1] for x in new_cigar if x[0]==0])

                prev_cigar_len = len(dict_copy[read_ID][read_num])-1
                prev_bp_aligned = 0
                parts = [x for x in dict_copy[read_ID][read_num].keys() if 'pt' in x]
                for part in parts:
                    if dict_copy[read_ID][read_num][part].cigar == 'match':
                        prev_bp_aligned += dict_copy[read_ID][read_num][part].ref_end-dict_copy[read_ID][read_num][part].ref_start

                # if read wasn't aligned before or CIGAR is shorter or cumulatively more aligned
                if len(new_cigar) < prev_cigar_len or new_bp_aligned > prev_bp_aligned+5:
                    idx = 1
                    nt_count = 0
                    nt_mapped = 0
                    for operation, length in read.cigar:
                        if operation in [0,1,4]:
                            if operation == 0: # matches
                                temp_dict['pt%i' % idx] = read_part(cigar='match',
                                                                    read_start=nt_count,
                                                                    read_end=nt_count+length,
                                                                    chrom=read.reference_name,
                                                                    strand=get_orientation(read.is_reverse),
                                                                    ref_start=read.pos+nt_mapped,
                                                                    ref_end=read.pos+nt_mapped+length,
                                                                    mapq=read.mapq)
                                nt_count += length
                                nt_mapped += length
                            elif operation == 1: # insertion
                                temp_dict['pt%i' % idx] = read_part(cigar='insertion',
                                                                    read_start=nt_count,
                                                                    read_end=nt_count+length)
                                nt_count += length
                            elif operation == 4: # soft-clipped
                                temp_dict['pt%i' % idx] = read_part(cigar='unmapped',
                                                                    read_start=nt_count,
                                                                    read_end=nt_count+length)
                                nt_count += length
                            idx += 1 # move to next index in dict if portion of read is clipped, aligned, or an insertion
                        elif operation == 2: # deletion
                            nt_mapped += length

                        # replace alignment in dict
                        dict_copy[read_ID][read_num] = temp_dict
                        
    return dict_copy


def align_clips(dict_in, ref_file, transgene_refname, fasta_name='temp.fa', sam_name='temp.sam', rDNA_refname='KY962518.1'):
    dict_copy = deepcopy(dict_in)

    seqs2write = []
    for matepair in dict_copy:
        for read in range(2):
            clips = [x for x in dict_copy[matepair][read].keys() if 'pt' in x]
            for clip in clips:
                end = dict_copy[matepair][read][clip].read_end
                start = dict_copy[matepair][read][clip].read_start
                if dict_copy[matepair][read][clip].cigar == 'unmapped' and (end-start) > 20:
                    seqs2write.append([matepair,str(read)+'.'+clip[-1],dict_copy[matepair][read]['seq'][start:end]])

    aligned_clips = get_alignments(seqs2write, ref_file, fasta_name, sam_name)

    for read in aligned_clips.fetch():
        read_ID = read.query_name[:-4]
        read_num = int(read.query_name[-3])
        clip_num = int(read.query_name[-1])
        total_clips = len([x for x in dict_copy[read_ID][read_num].keys() if 'pt' in x])
        temp_dict = {}

        temp_dict['seq'] = dict_copy[read_ID][read_num]['seq']
        if read.cigar and not read.is_supplementary:
            
            # determine whether ref was flipped
            strands = {}
            for x in dict_copy[read_ID][read_num].keys():
                if 'pt' in x and int(x[-1]) != clip_num and dict_copy[read_ID][read_num][x].cigar != 'unmapped':
                    strands[dict_copy[read_ID][read_num][x].chrom] = dict_copy[read_ID][read_num][x].strand
            flipped_ref = False
            if len(strands.values()) > 0: # False by default if length is zero
                flipped_ref = ('+' not in strands.values() and '-' in strands.values()) # if '-' in dict, but not '+', then true
                if '+' in strands.values() and '-' in strands.values():
                    if transgene_refname in strands:
                        flipped_ref = strands[transgene_refname] == '-'
                    elif rDNA_refname in strands:
                        flipped_ref = strands[rDNA_refname] == '-'
            strand = get_orientation(read.is_reverse ^ flipped_ref)

            # take all previous clips and add to temp dict
            if clip_num > 1:
                for x in range(1,clip_num):
                    temp_dict['pt%i' % x] = dict_copy[read_ID][read_num]['pt%i' % x]

            # if new portion is fully mapped, just swap out current entry
            if len(read.cigartuples) == 1:
                
                prev_entry = dict_copy[read_ID][read_num]['pt%i' % clip_num]
                temp_dict['pt%i' % clip_num] = read_part(cigar='match',
                                                         read_start=prev_entry.read_start,
                                                         read_end=prev_entry.read_end,
                                                         chrom=read.reference_name,
                                                         strand=strand,
                                                         ref_start=read.pos,
                                                         ref_end=read.pos+len(read.seq),
                                                         mapq=read.mapq)
                
                # take remaining clips and shift them
                if clip_num < total_clips:
                    for x in range(clip_num,total_clips):
                        temp_dict['pt%i' % (x+1)] = dict_copy[read_ID][read_num]['pt%i' % (x+1)]

            # if new portion has 2 or more parts, shift remaining portions down
            elif len(read.cigartuples) > 1:
                idx = clip_num
                nt_count = dict_copy[read_ID][read_num]['pt%i' % clip_num].read_start
                nt_mapped = 0
                if read.is_reverse:
                    for operation, length in read.cigar[::-1]:
                        if operation in [0,1,4]:
                            if operation == 0: # matches
                                temp_dict['pt%i' % idx] = read_part(cigar='match',
                                                                    read_start=nt_count,
                                                                    read_end=nt_count+length,
                                                                    chrom=read.reference_name,
                                                                    strand=strand,
                                                                    ref_start=read.pos+nt_mapped,
                                                                    ref_end=read.pos+nt_mapped+length,
                                                                    mapq=read.mapq)
                                nt_count += length
                                nt_mapped += length
                            elif operation == 1: # insertion
                                temp_dict['pt%i' % idx] = read_part(cigar='insertion',
                                                                    read_start=nt_count,
                                                                    read_end=nt_count+length)
                                nt_count += length
                            elif operation == 4: # soft-clipped
                                temp_dict['pt%i' % idx] = read_part(cigar='unmapped',
                                                                    read_start=nt_count,
                                                                    read_end=nt_count+length)
                                nt_count += length
                            idx += 1 # move to next index in dict if portion of read is clipped, aligned, or an insertion
                        elif operation == 2: # deletion
                            nt_mapped += length
                else:
                    for operation, length in read.cigar:
                        if operation in [0,1,4]:
                            if operation == 0: # matches
                                temp_dict['pt%i' % idx] = read_part(cigar='match',
                                                                    read_start=nt_count,
                                                                    read_end=nt_count+length,
                                                                    chrom=read.reference_name,
                                                                    strand=strand,
                                                                    ref_start=read.pos+nt_mapped,
                                                                    ref_end=read.pos+nt_mapped+length,
                                                                    mapq=read.mapq)
                                nt_count += length
                                nt_mapped += length
                            elif operation == 1: # insertion
                                temp_dict['pt%i' % idx] = read_part(cigar='insertion',
                                                                    read_start=nt_count,
                                                                    read_end=nt_count+length)
                                nt_count += length
                            elif operation == 4: # soft-clipped
                                temp_dict['pt%i' % idx] = read_part(cigar='unmapped',
                                                                    read_start=nt_count,
                                                                    read_end=nt_count+length)
                                nt_count += length
                            idx += 1 # move to next index in dict if portion of read is clipped, aligned, or an insertion
                        elif operation == 2: # deletion
                            nt_mapped += length

                # take remaining clips and shift them
                if clip_num < total_clips:
                    delta = idx - clip_num
                    for x in range(clip_num,total_clips):
                        temp_dict['pt%i' % (x+delta)] = dict_copy[read_ID][read_num]['pt%i' % (x+1)]

            dict_copy[read_ID][read_num] = temp_dict
                    
    return dict_copy


def filter_transgenemapping(dict_in, transgene_start, transgene_end, transgene_refname, buffer=5):
    df_out = {}
    for matepair in dict_in:
        matepair_maps = False
        for read in range(2):
            parts = [x for x in dict_in[matepair][read].keys() if 'pt' in x]
            for part in parts:
                matepair_maps = matepair_maps or \
                    (dict_in[matepair][read][part].chrom == transgene_refname and \
                     (transgene_end > dict_in[matepair][read][part].ref_start > transgene_start or \
                      transgene_end > dict_in[matepair][read][part].ref_end > transgene_start+buffer))
        if matepair_maps:
            df_out[matepair] = dict_in[matepair]
            
    return df_out


# try again to map the unaligned clipped portions of reads
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

def remap_sensitive(dict_in, template_sequence, rDNA_sequence, transgene_refname,
                    TS_context='TGACTCTCTTAAGGTAGCCAAA', TS_offset=11665, rDNA_refname='KY962518.1'):
    dict_copy = deepcopy(dict_in)
    
    for matepair in dict_copy:
        for read in range(2):
            parts = [x for x in dict_copy[matepair][read].keys() if 'pt' in x]
            for part in parts:
                prev_entry = dict_copy[matepair][read][part]
                seq_len = prev_entry.get_length()
                seq = dict_copy[matepair][read]['seq'][prev_entry.read_start:prev_entry.read_end]
                
                # first, check if it's a substring of known target site
                if prev_entry.cigar == 'unmapped' and seq_len > 5:
                    
                    # determine whether ref was flipped
                    strands = {}
                    for x in parts:
                        if x != part and dict_copy[matepair][read][x].cigar != 'unmapped':
                            strands[dict_copy[matepair][read][x].chrom] = dict_copy[matepair][read][x].strand
                    flipped_ref = False
                    if len(strands.values()) > 0: # False by default if length is zero
                        flipped_ref = ('+' not in strands.values() and '-' in strands.values()) # if '-' in dict, but not '+', then true
                        if '+' in strands.values() and '-' in strands.values():
                            if transgene_refname in strands:
                                flipped_ref = strands[transgene_refname] == '-'
                            elif rDNA_refname in strands:
                                flipped_ref = strands[rDNA_refname] == '-'
                    
                    if seq in TS_context:
                        strand = get_orientation(not flipped_ref)
                        dict_copy[matepair][read][part] = read_part(cigar='match',
                                                                    read_start=prev_entry.read_start,
                                                                    read_end=prev_entry.read_end,
                                                                    chrom=rDNA_refname,
                                                                    strand=strand,
                                                                    ref_start=TS_context.find(seq)+TS_offset,
                                                                    ref_end=TS_context.find(seq)+TS_offset+seq_len)
                    elif rc(seq) in TS_context:
                        strand = get_orientation(flipped_ref)
                        dict_copy[matepair][read][part] = read_part(cigar='match',
                                                                    read_start=prev_entry.read_start,
                                                                    read_end=prev_entry.read_end,
                                                                    chrom=rDNA_refname,
                                                                    strand=strand,
                                                                    ref_start=TS_context.find(rc(seq))+TS_offset,
                                                                    ref_end=TS_context.find(rc(seq))+TS_offset+seq_len)
                        
                # then check 45S and transgene again
                elif dict_copy[matepair][read][part].cigar == 'unmapped' and seq_len > 12:
                    dist = int(seq_len/10)
                    
                    # next align to 45S or transgene sequence
                    for sequence, name in [(rDNA_sequence, rDNA_refname), (template_flanks, transgene_refname)]:
                        match = get_bestmatch(find_near_matches(seq, rDNA_sequence, max_l_dist=dist))
                        if match != 0:
                            strand = get_orientation(not flipped_ref)
                            dict_copy[matepair][read][part] = read_part(cigar='match',
                                                                        read_start=prev_entry.read_start,
                                                                        read_end=prev_entry.read_end,
                                                                        chrom=rDNA_refname,
                                                                        strand=strand,
                                                                        ref_start=match.start,
                                                                        ref_end=match.end)
                            break
                        else:
                            match = get_bestmatch(find_near_matches(rc(seq), rDNA_sequence, max_l_dist=dist))
                            if match != 0:
                                strand = get_orientation(flipped_ref)
                                dict_copy[matepair][read][part] = read_part(cigar='match',
                                                                        read_start=prev_entry.read_start,
                                                                        read_end=prev_entry.read_end,
                                                                        chrom=rDNA_refname,
                                                                        strand=strand,
                                                                        ref_start=match.start,
                                                                        ref_end=match.end)
                                break                    
    
    return dict_copy

def is_rDNA(chrom, coord):
    chr13 = chrom == 'chr13' and 5770548 < coord < 9348041
    chr14 = chrom == 'chr14' and 2099537 < coord < 2817811
    chr15 = chrom == 'chr15' and 2506442 < coord < 4707485
    chr21 = chrom == 'chr21' and 3108298 < coord < 5612715
    chr22 = chrom == 'chr22' and 4793794 < coord < 5720650
    is_rDNA = chr13 or chr14 or chr15 or chr21 or chr22
    return is_rDNA

def realign_rDNA_reads(dict_in, transgene_refname, tg_start, tg_end, rDNA_sequence, rDNA_refname='KY962518.1'):
    
    dict_copy = deepcopy(dict_in)
    
    for matepair in dict_copy:
        for read in range(2):
            parts = [x for x in dict_copy[matepair][read].keys() if 'pt' in x]
            for part in parts:
                entry = dict_copy[matepair][read][part]
                if is_rDNA(entry.chrom, entry.ref_end) or (entry.chrom==transgene_refname and (entry.ref_start>tg_end-2 or entry.ref_end < tg_start+2)):
                    seq_len = entry.get_length()
                    seq = dict_copy[matepair][read]['seq'][entry.read_start:entry.read_end]
                    dist = int(seq_len/10)
                    match = get_bestmatch(find_near_matches(seq, rDNA_sequence, max_l_dist=dist))
                    if match != 0:
                        dict_copy[matepair][read][part] = read_part(cigar='match',
                                                                    read_start=entry.read_start,
                                                                    read_end=entry.read_end,
                                                                    chrom=rDNA_refname,
                                                                    strand=entry.strand,
                                                                    ref_start=match.start,
                                                                    ref_end=match.end)
                        break
                    else:
                        match = get_bestmatch(find_near_matches(rc(seq), rDNA_sequence, max_l_dist=dist))
                        if match != 0:
                            dict_copy[matepair][read][part] = read_part(cigar='match',
                                                                        read_start=entry.read_start,
                                                                        read_end=entry.read_end,
                                                                        chrom=rDNA_refname,
                                                                        strand=entry.strand,
                                                                        ref_start=match.start,
                                                                        ref_end=match.end)
                            break                    
    
    return dict_copy

def discard_plasmid_reads(dict_in, ref_file, fasta_name='temp.fa', sam_name='temp.sam'):
    dict_copy = deepcopy(dict_in)
    
    seqs2write = []
    for matepair in dict_copy:
        for read in range(2):
            # if not fully aligned in one part
            if not (dict_copy[matepair][read]['pt1'].cigar == 'match' and 'pt2' not in dict_copy[matepair][read]):
                seqs2write.append([matepair,read,dict_copy[matepair][read]['seq']])

    plasmid_alignments = get_alignments(seqs2write, ref_file, fasta_name, sam_name)

    for read in plasmid_alignments.fetch():
        read_ID = read.query_name[:-2]
        read_num = int(read.query_name[-1])
        temp_dict = {}
        if len(read.cigar) == 1 and read_ID in dict_copy:
            del dict_copy[read_ID]
                        
    return dict_copy

def discard_contaminating_reads(dict_in, ref_file, fasta_name='temp.fa', sam_name='temp.sam'):
    dict_copy = deepcopy(dict_in)
    seqs2write = []
    for matepair in dict_in:
        for read in range(2):
            # if not fully aligned in one part
            if not (dict_copy[matepair][read]['pt1'].cigar == 'match' and 'pt2' not in dict_copy[matepair][read]):
                seqs2write.append([matepair,read,dict_copy[matepair][read]['seq']])

    blacklist_alignments = get_alignments(seqs2write, ref_file, fasta_name, sam_name)

    for read in blacklist_alignments.fetch():
        read_ID = read.query_name[:-2]
        if read.cigar and read_ID in dict_copy and read.alen>60:
            del dict_copy[read_ID]
    
    return dict_copy
