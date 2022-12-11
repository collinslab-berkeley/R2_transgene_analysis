import numpy as np
import pandas as pd
import pysam
from Bio.Seq import Seq
from Bio import SeqIO
import os, sys

# GLOBAL VARIABLES
COLUMNS = ['read_ID', 'read1_orientation', 'read2_orientation', 'read1_seq', 'read2_seq',
           'read1_chr1', 'read1_start1', 'read1_end1', 'read1_loc1', 'read1_mapq1', 'read1_clipstatus1',
           'read1_chr2', 'read1_start2', 'read1_end2', 'read1_loc2', 'read1_mapq2', 'read1_clipstatus2',
           'read2_chr1', 'read2_start1', 'read2_end1', 'read2_loc1', 'read2_mapq1', 'read2_clipstatus1',
           'read2_chr2', 'read2_start2', 'read2_end2', 'read2_loc2', 'read2_mapq2', 'read2_clipstatus2']


read_11_cols = ['read1_chr1', 'read1_start1', 'read1_end1', 'read1_loc1', 'read1_mapq1', 'read1_clipstatus1']
read_12_cols = ['read1_chr2', 'read1_start2', 'read1_end2', 'read1_loc2', 'read1_mapq2', 'read1_clipstatus2']
read_21_cols = ['read2_chr1', 'read2_start1', 'read2_end1', 'read2_loc1', 'read2_mapq1', 'read2_clipstatus1']
read_22_cols = ['read2_chr2', 'read2_start2', 'read2_end2', 'read2_loc2', 'read2_mapq2', 'read2_clipstatus2']

col_dict = {'read1.1': read_11_cols, 'read1.2': read_12_cols, 'read2.1': read_21_cols, 'read2.2': read_22_cols}

rDNA_chroms = ['chr13','chr14','chr15','chr21','chr22']

def classify_read(cigar):
    if cigar == None:
        return 'unmapped'
    elif len(cigar) == 1:
        return 'fully_mapped'
    elif cigar[0][0] == 4:
        if cigar[-1][0] == 4:
            return 'middle_mapped'
        else:
            return '5clipped'
    else:
        return '3clipped'
    
def get_orientation(reverse_bool):
    if reverse_bool:
        return 'R'
    else:
        return 'F'

def update_data(read, clip_status, entry):
    if read.is_read1:
        if clip_status == 'fully_mapped' and not (read.reference_name in rDNA_chroms and read.mapq < 10):
            entry = [get_orientation(read.is_reverse)] + [entry[1]] + [read.seq] + [entry[3]] + \
                    [read.reference_name, read.pos, read.aend, 0, read.mapq, clip_status]+entry[-18:]
        elif clip_status == '3clipped':
            if read.reference_name in rDNA_chroms and read.mapq < 10:
                pass
            elif entry[4] == -1:
                entry = [get_orientation(read.is_reverse)] + [entry[1]] + [read.seq] + [entry[3]] + \
                        [read.reference_name, read.pos, read.aend, 0, read.mapq, clip_status]+entry[-18:]
            elif read.seq == str(Seq(entry[2]).reverse_complement()):
                loc = read.cigartuples[-1][1]
                entry = entry[:10] + [read.reference_name, read.pos, read.aend, loc, read.mapq, clip_status]+entry[-12:]
            # else:
            #     print('some other issue')
            #     print(entry)
            #     print(read)
            #     print(' ')
        elif clip_status == '5clipped':
            loc = read.cigartuples[0][1]
            if read.reference_name in rDNA_chroms and read.mapq < 10:
                pass
            elif entry[10] == -1:
                entry = [get_orientation(read.is_reverse)] + [entry[1]] + [read.seq] + entry[3:10] + \
                    [read.reference_name, read.pos, read.aend, loc, read.mapq, clip_status]+entry[-12:]
            elif read.seq == str(Seq(entry[2]).reverse_complement()):
                entry = entry[:4] + [read.reference_name, read.pos, read.aend, 0, read.mapq, clip_status]+entry[-18:]
            # else:
            #     print('some other issue')
            #     print(entry)
            #     print(read)
            #     print(' ')
    else:
        if clip_status == 'fully_mapped' and not (read.reference_name in rDNA_chroms and read.mapq < 10):
            entry = [entry[0]] + [get_orientation(read.is_reverse)] + [entry[2]] + [read.seq] + entry[4:16] + \
                    [read.reference_name, read.pos, read.aend, 0, read.mapq, clip_status]+entry[-6:]
        elif clip_status == '3clipped':
            if read.reference_name in rDNA_chroms and read.mapq < 10:
                pass
            elif entry[16] == -1:
                entry = [entry[0]] + [get_orientation(read.is_reverse)] + [entry[2]] + [read.seq] + entry[4:16] + \
                    [read.reference_name, read.pos, read.aend, 0, read.mapq, clip_status]+entry[-6:]
            elif read.seq == str(Seq(entry[3]).reverse_complement()):
                loc = read.cigartuples[-1][1]
                entry = entry[:22] + [read.reference_name, read.pos, read.aend, loc, read.mapq, clip_status]
            # else:
            #     print('some other issue')
            #     print(entry)
            #     print(read)
            #     print(' ')
        elif clip_status == '5clipped':
            loc = read.cigartuples[0][1]
            if read.reference_name in rDNA_chroms and read.mapq < 10:
                pass
            elif entry[22] == -1:
                entry = [entry[0]] + [get_orientation(read.is_reverse)] + [entry[2]] + [read.seq] + entry[4:22] + \
                    [read.reference_name, read.pos, read.aend, loc, read.mapq, clip_status]
            elif read.seq == str(Seq(entry[3]).reverse_complement()):
                entry = entry[:16] + [read.reference_name, read.pos, read.aend, 0, read.mapq, clip_status]+entry[-6:]
            # else:
            #     print('some other issue')
            #     print(entry)
            #     print(read)
            #     print(' ')
                
    return entry

def add_unmapped_read(read, entry):
    if read.is_read1:
        entry = entry[:2]+[read.seq]+entry[3:]
    else:
        entry = entry[:3]+[read.seq]+entry[4:]
    return entry

def realign_reads(df, transgene_refname, fasta_name, sam_name, ref_file):
    temp_df = df.copy()
    seqs2write = []
    for i, row in temp_df.iterrows():
        if row['read1_orientation']==-1 and row['read1_seq']!=-1:
            seqs2write.append([row['read_ID'],1,row['read1_seq']])
        elif row['read1_clipstatus1']=='3clipped' and row['read1_chr1']==transgene_refname:
            seqs2write.append([row['read_ID'],1,row['read1_seq']])
        elif row['read1_clipstatus2']=='5clipped' and row['read1_chr2']==transgene_refname:
            seqs2write.append([row['read_ID'],1,row['read1_seq']])
        if row['read2_orientation']==-1 and row['read2_seq']!=-1:
            seqs2write.append([row['read_ID'],2,row['read2_seq']])
        elif row['read2_clipstatus1']=='3clipped' and row['read2_chr1']==transgene_refname:
            seqs2write.append([row['read_ID'],2,row['read2_seq']])
        elif row['read2_clipstatus2']=='5clipped' and row['read2_chr2']==transgene_refname:
            seqs2write.append([row['read_ID'],2,row['read2_seq']])

    with open(fasta_name, 'w') as fa_out:
        for seq in seqs2write:
            fa_out.write('>%s_%i\n' % (seq[0],seq[1]))
            fa_out.write(seq[2]+'\n')
    
    os.system('bwa mem %s %s > %s' % (ref_file, fasta_name, sam_name))
    
    realigned_reads = pysam.AlignmentFile(sam_name, 'r')

    for read in realigned_reads.fetch():
        clip_status = classify_read(read.cigartuples)
        read_ID = read.query_name[:-2]
        read_num = int(read.query_name[-1])
        loc = temp_df[temp_df['read_ID']==read_ID].index[0]
        if clip_status == 'fully_mapped':
            if read_num == 1:
                temp_df.loc[loc,read_11_cols] = [read.reference_name,read.pos,read.aend,0,read.mapq,clip_status]
                temp_df.loc[loc,read_12_cols] = -1
            else:
                temp_df.loc[loc,read_21_cols] = [read.reference_name,read.pos,read.aend,0,read.mapq,clip_status]
                temp_df.loc[loc,read_22_cols] = -1
            if read.is_reverse:
                temp_df.loc[loc,'read%i_orientation' % read_num] = 'R'
                temp_df.loc[loc,'read%i_seq' % read_num] = str(Seq(read.seq).reverse_complement())
            else:
                temp_df.loc[loc,'read%i_orientation' % read_num] = 'F'
                temp_df.loc[loc,'read%i_seq' % read_num] = read.seq
        elif clip_status == '5clipped' and (read.alen-5) > (temp_df.loc[loc,'read%i_end2' % read_num]-temp_df.loc[loc,'read%i_start2' % read_num]):
            start = read.cigartuples[0][1]
            if read_num == 1 :
                temp_df.loc[loc,read_12_cols] = [read.reference_name,read.pos,read.aend,start,read.mapq,clip_status]
            else:
                temp_df.loc[loc,read_22_cols] = [read.reference_name,read.pos,read.aend,start,read.mapq,clip_status]
            if temp_df.loc[loc,'read%i_orientation' % read_num] == -1:
                if read.is_reverse:
                    temp_df.loc[loc,'read%i_orientation' % read_num] = 'R'
                    temp_df.loc[loc,'read%i_seq' % read_num] = str(Seq(read.seq).reverse_complement())
                else:
                    temp_df.loc[loc,'read%i_orientation' % read_num] = 'F'
                    temp_df.loc[loc,'read%i_seq' % read_num] = read.seq
        elif clip_status == '3clipped' and (read.alen-5) > (temp_df.loc[loc,'read%i_end1' % read_num]-temp_df.loc[loc,'read%i_start1' % read_num]):
            if read_num == 1 :
                temp_df.loc[loc,read_11_cols] = [read.reference_name,read.pos,read.aend,0,read.mapq,clip_status]
            else:
                temp_df.loc[loc,read_21_cols] = [read.reference_name,read.pos,read.aend,0,read.mapq,clip_status]
            if temp_df.loc[loc,'read%i_orientation' % read_num] == -1:
                if read.is_reverse:
                    temp_df.loc[loc,'read%i_orientation' % read_num] = 'R'
                    temp_df.loc[loc,'read%i_seq' % read_num] = str(Seq(read.seq).reverse_complement())
                else:
                    temp_df.loc[loc,'read%i_orientation' % read_num] = 'F'
                    temp_df.loc[loc,'read%i_seq' % read_num] = read.seq

    return temp_df

def has_alignment_gap(df_row, read):
    return (df_row['read%i_end1' % read] - df_row['read%i_start1' % read] + \
            df_row['read%i_end2' % read] - df_row['read%i_start2' % read]) < (len(df_row['read%i_seq' % read]) - 5)

def align_clips(df, transgene_refname, fasta_name, sam_name, ref_file):
    temp_df = df.copy()
    seqs2write = []
    for i, row in temp_df.iterrows():
        # read 1
        start = int(row['read1_end1']-row['read1_start1'])
        end = int(row['read1_loc2'])
        if row['read1_clipstatus1']=='3clipped' and row['read1_clipstatus2']!='5clipped':
            seqs2write.append([row['read_ID'],1,2,row['read1_seq'][start:]])
        elif row['read1_clipstatus2']=='5clipped' and row['read1_clipstatus1']!='3clipped':
            seqs2write.append([row['read_ID'],1,1,row['read1_seq'][:end]])
        elif has_alignment_gap(row, 1):
            seqs2write.append([row['read_ID'],1,1,row['read1_seq'][:end]])
            seqs2write.append([row['read_ID'],1,2,row['read1_seq'][start:]])
        # read 2
        start = int(row['read2_end1']-row['read2_start1'])
        end = int(row['read2_loc2'])
        if row['read2_clipstatus1']=='3clipped' and row['read2_clipstatus2']!='5clipped':
            seqs2write.append([row['read_ID'],2,2,row['read2_seq'][start:]])
        elif row['read2_clipstatus2']=='5clipped' and row['read2_clipstatus1']!='3clipped':
            seqs2write.append([row['read_ID'],2,1,row['read2_seq'][:end]])
        elif has_alignment_gap(row, 2):
            seqs2write.append([row['read_ID'],2,1,row['read2_seq'][:end]])
            seqs2write.append([row['read_ID'],2,2,row['read2_seq'][start:]])
    
    seqs2write = [x for x in seqs2write if len(x[3])>29]
            
    with open(fasta_name, 'w') as fa_out:
        for seq in seqs2write:
            fa_out.write('>%s_%i.%i\n' % (seq[0],seq[1],seq[2]))
            fa_out.write(seq[3]+'\n')
    
    os.system('bwa mem %s %s > %s' % (ref_file, fasta_name, sam_name))
    
    realigned_reads = pysam.AlignmentFile(sam_name, 'r')

    for read in realigned_reads.fetch():
        clip_status = classify_read(read.cigartuples)
        read_ID = read.query_name[:-4]
        read_num = int(read.query_name[-3])
        read_end = int(read.query_name[-1])
        loc = temp_df[temp_df['read_ID']==read_ID].index[0]
        
        # if clip maps more than current read, then replace
        if clip_status == 'fully_mapped' or (read.alen and (read.alen - 5) > (temp_df.loc[loc,'read%i_end%i' % (read_num, read_end)] - \
                                                                             temp_df.loc[loc,'read%i_start%i' % (read_num, read_end)])):
            if read_end == 1:
                if read.is_reverse:
                    clip = '5clipped'
                    start = temp_df.loc[loc,'read%i_seq' % read_num].find(str(Seq(read.query_sequence).reverse_complement()))
                else:
                    clip = '3clipped'
                    start = temp_df.loc[loc,'read%i_seq' % read_num].find(read.query_sequence)
                
            elif read_end == 2:
                if read.is_reverse:
                    clip = '3clipped'
                    start = temp_df.loc[loc,'read%i_seq' % read_num].find(str(Seq(read.query_sequence).reverse_complement()))
                else:
                    clip = '5clipped'
                    start = temp_df.loc[loc,'read%i_seq' % read_num].find(read.query_sequence)
            temp_df.loc[loc,col_dict['read%i.%i' % (read_num,read_end)]] = [read.reference_name,read.pos,read.aend,start,read.mapq,clip]
    return temp_df

def filter_transgenemapping(df, flank_len, transgene_len, transgene_refname, offset=5):
    temp_df = df.copy()
    transgene_start = flank_len
    transgene_end = flank_len+transgene_len

    maps_to_transgene = []
    for i, row in temp_df.iterrows():
        # print(row)
        row_maps = False
        for x,y in zip([1,1,2,2],[1,2,1,2]):
            row_maps = row_maps or \
                        (row['read%i_chr%i' % (x,y)] == transgene_refname and \
                        ((transgene_end-offset) > row['read%i_start%i' % (x,y)] > (transgene_start+offset) or \
                         (transgene_end-offset) > row['read%i_end%i' % (x,y)] > (transgene_start+offset)))
        maps_to_transgene.append(row_maps)
        
    return temp_df[maps_to_transgene].reset_index(drop=True)


def main():
    transgene_ref, input_sam, output_csv = sys.argv[1:]
    flank_len = 840
    transgene_len = len(list(SeqIO.parse(transgene_ref, format='fasta'))[0].seq) - 2*flank_len
    transgene_refname = list(SeqIO.parse(transgene_ref, format='fasta'))[0].name

    # read sam file and convert to dataframe representation
    transgenereads = pysam.AlignmentFile(input_sam, 'r')
    read_dict = {}
    for read in transgenereads.fetch():
        clip_status = classify_read(read.cigartuples)
        if read.query_name and not read.is_supplementary:
            if clip_status in ['fully_mapped', '5clipped', '3clipped'] and read.alen > 30:
                if read.query_name in read_dict:
                    read_dict[read.query_name] = update_data(read, clip_status, read_dict[read.query_name])
                else:
                    read_dict[read.query_name] = update_data(read, clip_status, [-1]*28)
            elif read.query_name in read_dict:
                read_dict[read.query_name] = add_unmapped_read(read, read_dict[read.query_name])
            else:
                read_dict[read.query_name] = add_unmapped_read(read, [-1]*28)
                
    # create DF
    read_df = pd.DataFrame.from_records(read_dict).T.reset_index()
    read_df.columns = COLUMNS

    # filter to make sure at least one read crosses transgene
    read_df = filter_transgenemapping(read_df, flank_len, transgene_len, transgene_refname)

    # test better alignment elsewhere (in 45S rDNA and across genome) 
    read_df = realign_reads(read_df, transgene_refname, 'temp.fa', 'temp.sam', 'references/rDNAscaffold_45S.fa')
    read_df = realign_reads(read_df, transgene_refname, 'temp.fa', 'temp.sam', 'references/chm13v2.0.fa.gz')
    read_df = filter_transgenemapping(read_df, flank_len, transgene_len, transgene_refname)

    # filter by proper read pairs (FR or RF orientation; both pairs mapped)
    read_df['mate_orientation'] = read_df['read1_orientation'].map(str) + ',' + read_df['read2_orientation'].map(str)
    read_df = read_df[read_df['mate_orientation'].isin(['R,F','F,R'])]

    # align clipped portions of reads
    read_df = align_clips(read_df, transgene_refname, 'temp.fa', 'temp.sam', transgene_ref)
    read_df = align_clips(read_df, transgene_refname, 'temp.fa', 'temp.sam', 'references/rDNAscaffold_45S.fa')
    read_df = align_clips(read_df, transgene_refname, 'temp.fa', 'temp.sam', 'references/chm13v2.0.fa.gz')

    # save to file for downstream analysis
    read_df = read_df.reset_index(drop=True)
    read_df.to_csv(output_csv)


if __name__ == '__main__':
    main()