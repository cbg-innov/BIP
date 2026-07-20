#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 16:02:45 2025
Includes "second-tier check": for sequences that still contain indels after homopolymer-based
corrections, we attempt to correct indels just based on alignment to the reference sequence. 
@author: robinfloyd
"""
import sys
import csv
import os
import subprocess
from contextlib import chdir
from Bio import Align, SeqIO
from Bio.Seq import Seq

def match_by_vsearch(qry_file,reference_set):
    hit_seqs_dict = {}
    vsearch_dict = {}
    output_file = 'vsearch_hit.tsv'

    vsearch_command = [ 'vsearch',
                       '--usearch_global',
                       qry_file,
                       '--db',
                       reference_set,
                       '--blast6out',
                       output_file,
                       '--id', '0.7',
                       '--maxaccepts', '1',
                       '--maxhits', '1'
                       ]

    subprocess.run(vsearch_command)

    with open(output_file, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            query_name, hit_name = line[:2]
            vsearch_dict[query_name] = hit_name
    
    top_hits = set(vsearch_dict.values())
    
    for record in SeqIO.parse(reference_set, 'fasta'):
        if record.id in top_hits:
            hit_seqs_dict[record.id] = str(record.seq.upper())
        
    return vsearch_dict, hit_seqs_dict

def end_correction(ref_aln, qry_aln):
    
    ref_aln_endfix = ref_aln
    qry_aln_endfix = qry_aln
    endfixed = False
    startgapcounter = 0
    endgapcounter = 0
    
    # If reference sequence starts or ends with gaps, query sequence is too long.
    # Count how many consecutive gaps, trim those characters from both sequences at either or both ends as needed.
    if ref_aln_endfix.startswith('-'): 
        counter = 0   
        while ref_aln_endfix[counter] == '-':
            counter += 1
        ref_aln_endfix = ref_aln_endfix[counter:]
        qry_aln_endfix = qry_aln_endfix[counter:]
        endfixed = True
        # Also need to keep track of how many characters we removed from the query sequence,
        # so we can output them later correctly aligned with the original sequences 
        startgapcounter = counter
    
    if ref_aln_endfix.endswith('-'):
        counter = (len(ref_aln_endfix)-1) # because indexing starts at 0, the index of the last position is the length of the string -1
        # Also need to keep track of how many characters we removed from the query sequence,
        # so we can output them later correctly aligned with the original sequences 
        endgapcounter = 0
        while ref_aln_endfix[counter] == '-':
            counter -= 1
            endgapcounter += 1
        ref_aln_endfix = ref_aln_endfix[:counter+1]
        qry_aln_endfix = qry_aln_endfix[:counter+1]
        endfixed = True
    
    # If query sequence starts or ends with gaps, replace them with Ns.
    if qry_aln_endfix.startswith('-'): 
        counter = 0   
        while qry_aln_endfix[counter] == '-':
            counter += 1
        qry_aln_endfix = ('n' * counter) + qry_aln_endfix[counter:]
        endfixed = True
    
    if qry_aln_endfix.endswith('-'):
        counter = (len(qry_aln_endfix)-1) # because indexing starts at 0, the index of the last position is the length of the string -1
        while qry_aln_endfix[counter] == '-':
            counter -= 1
        qry_aln_endfix = qry_aln_endfix[:counter+1] + ('n' * (len(qry_aln_endfix) - counter - 1))
        endfixed = True

    return ref_aln_endfix, qry_aln_endfix, endfixed, startgapcounter, endgapcounter

def homopolymer_length(seq, index):
    """
    Given a sequence and a position index (where a gap occurs),
    return the length of the homopolymer run around that index.
    The index refers to the position in the *aligned* reference.
    """
    if index < 0 or index >= len(seq):
        return '','','','',''

    # Determine the character
    base = seq[index]
    left_index = index - 1

    if left_index >= 0:
        while left_index >= 0 and seq[left_index] == "-":
            left_index -= 1
        base_left = seq[left_index]
    else:
        base_left = ''
    
    last_gap_position = index
    right_index = index + 1

    if right_index < len(seq):
        while right_index < len(seq) and seq[right_index] == "-":
            last_gap_position = right_index
            right_index += 1
        base_right = seq[right_index]
    else:
        base_right = ''

    # Count to left
    left_count = 0
    l = left_index
    while l >= 0 and (seq[l] == base_left or seq[l] == '-' or seq[l] == 'N'):
        left_start_index = l
        if seq[l] != '-':
            left_count += 1
        l -= 1
        
    homopolymer_left = seq[(l+1):index].replace('-','')
    
    # Count to right
    right_count = 0
    r = right_index
    while r < len(seq) and (seq[r] == base_right or seq[r] == '-' or seq[r] == 'N'):
        right_end_index = r
        if seq[r] != '-':
            right_count += 1
        r += 1
    
    homopolymer_right = seq[right_index:r].replace('-','')
   
    if base == '-': # Possibility 1: position is a gap
    
        # Possibility 1a: position is gap, letters either side the same as each other
        if base_left == base_right:
            homopolymer = homopolymer_left + homopolymer_right
            placement = 'part'
            # left and right indices unchanged

        # Possibility 1b: position is gap, letters either side different from each other
    
        else:            
            if len(homopolymer_left) >= len(homopolymer_right):
                homopolymer = homopolymer_left
                right_end_index = index - 1
                # left index unchanged

            elif len(homopolymer_left) < len(homopolymer_right):
                homopolymer = homopolymer_right
                left_start_index = last_gap_position + 1  
                # right index unchanged
            placement = 'part'
            

    elif base == 'N': # Possibility 2: position is N, i.e. an unspecified base. Functionally this can be treated as part of a homopolymer.
        
        # Possibility 2a: position is N, letters either side the same as each other.
        # Add together with the N in the middle.
        if (base_left == base_right) or (base_left == 'N') or (base_right == 'N'):
            homopolymer = homopolymer_left + base + homopolymer_right
            placement = 'part'
            # left and right indices unchanged
        
        # Possibility 2b: position is N, letters on left and right are different from each other.
        # Add the N to whichever is longer. If they are the same length, arbitrarily add to the left one.
        else: 
            if len(homopolymer_left) >= len(homopolymer_right):
                homopolymer = homopolymer_left + base
                right_end_index = index
                # left index unchanged

            elif len(homopolymer_left) < len(homopolymer_right):
                homopolymer = base + homopolymer_right
                left_start_index = index  
                # right index unchanged

            placement = 'part'

    else: # Possibility 3: position is a specific base (i.e. letter other than N)
        
        # Possibility 3a: position is letter, letters either side the same as the position and each other
        if base == base_left == base_right:
            homopolymer = homopolymer_left + base + homopolymer_right
            placement = 'part'
            # left and right indices unchanged

       # Possibility 3b: position is letter, letter on left is the same, letter on right different
        elif base == base_left != base_right:
            homopolymer = homopolymer_left + base          
            placement = 'part'
            right_end_index = index - 1
            # left index unchanged
            
        # Possibility 3c: position is letter, letter on right is the same, letter on left different
        elif base == base_right != base_left:
            homopolymer = base + homopolymer_right
            placement = 'part'
            left_start_index = last_gap_position + 1  
            # right index unchanged
        
        # Possibility 3d: position is letter, letters on both left and right are different from it, irrelevant whether L & R same as each other
        # Take whichever of left or right is longer. If equal, arbitrarily take the left one.
        elif (base != base_right) & (base != base_left):
            if len(homopolymer_left) >= len(homopolymer_right):
                homopolymer = homopolymer_left
                right_end_index = index - 1
                # left index unchanged
                placement = 'adjacent-L'
            else:
                homopolymer = homopolymer_right
                placement = 'adjacent-R'
                left_start_index = last_gap_position + 1  
                # right index unchanged            
            
    return base, placement, homopolymer, left_start_index, right_end_index

def second_tier_check(seq_name, qry_aln_edited, indels_not_corrected_list):
    for run in indels_not_corrected_list:
        query_or_reference = run.pop(0)
        # Take the first letter of the list to tell us if this is a gap in the query or the reference.
        # This also removes the letter from the current gap run list.
                
        if len(run) <= 2: # only try to correct 1 or 2 consecutive gaps
            if query_or_reference == 'r':
                # Gap is in the reference sequence.
                # Delete the extra letters (replace with gaps)
                qry_aln_edited = qry_aln_edited[:run[0]] + ('-' * len(run)) + qry_aln_edited[(run[-1]+1):]
                print('Deleted bases from sequence', seq_name,'at position(s)', run, '(non-homopolymer)')

            elif query_or_reference == 'q':
                # Gap is in the query sequence.
                # Replace these gaps with Ns, and leave the reference alone.
                qry_aln_edited = qry_aln_edited[:run[0]] + ('n' * len(run)) + qry_aln_edited[(run[-1]+1):]
                print('Added Ns to sequence', seq_name,'at position(s)', run, '(non-homopolymer)')               

        return(qry_aln_edited)

def hmm_check(input_sequence, trans_table):
    sequence_nogaps = input_sequence.replace('-','')
    length_issue = False
    stop_codon = None

    # Python's inbuilt translator will throw an error if the sequence length is not a multiple of 3.
    # Next part pads the end of the sequence with Ns to avoid this
    while (len(sequence_nogaps)-1)%3 != 0:
        length_issue = True
        sequence_nogaps = sequence_nogaps + 'N'
    
    # Translate to amino acids, starting from the 2nd position in the sequence (aka index 1 in python)            
    protein_seq = str(Seq(sequence_nogaps[1:]).translate(table=trans_table))
    
    if '*' in protein_seq:
        stop_codon = True
    else:
        stop_codon = False

    return stop_codon, length_issue

# Main script starts here

aligner = Align.PairwiseAligner() # Alignment function from Biopython. Parameters set below
aligner.mode = 'global'

# # Core substitution scoring
aligner.match_score = 2
aligner.mismatch_score = -4

# # Internal gaps: strongly discouraged
aligner.open_internal_gap_score = -20
aligner.extend_internal_gap_score = -5

# # END GAPS — asymmetric on purpose
# # Allow trimming of QUERY only
aligner.open_end_deletion_score = 0
aligner.extend_end_deletion_score = 0

# # Do NOT allow free padding of query with Ns
aligner.open_end_insertion_score = -10
aligner.extend_end_insertion_score = -2

qry_file = sys.argv[1]
reference_set = sys.argv[2]

qry_seq_list = list(SeqIO.parse(qry_file, 'fasta'))
output_dict = {}
uncertain_dict = {}
problem_dict = {}

alignment_directory = "alignment_dir"
os.makedirs(alignment_directory, exist_ok=True)

vsearch_dict, hit_seqs_dict = match_by_vsearch(qry_file,reference_set)

output_no_change_count = 0
excluded_for_length_count = 0
excluded_for_nonmatch_count = 0
endfixed_count = 0
no_changes_except_endfix_count = 0
autocorrect_count = 0
non_hp_autocorrect_count = 0
indel_not_corrected_count = 0
stop_count = 0

for qry_record in qry_seq_list:
    
    seq_name = qry_record.id
    qry_seq = str(qry_record.seq.upper())
    qry_seq_length = len(qry_seq)
    
    endfixed = None
    length_issue = None
    no_issues = None
    stop_codon = None
    sequence_changed = None
    indels_not_corrected_list = []
    hmm_check_needed = True
    autocorrect_bool_list = []
    qry_aln_second_tier = None
    temp_name = None
    non_hp_edit = False
    startgapcounter = 0
    endgapcounter = 0

    if not seq_name in vsearch_dict.keys():
        # if the sequence found no match, exclude
        print('Sequence', seq_name, '- no match found.')
        edited_name = seq_name + '|no_match'
        problem_dict[edited_name] = qry_seq
        excluded_for_nonmatch_count += 1
        continue
    
    # Get the details of the top hit   
    hit = vsearch_dict[seq_name]
    ref_seq = hit_seqs_dict[hit]
    # Get codon translation table as the characters after the last '|' in the fasta header of the topt hit;
    # This is assumed to be the correct translation table for this query sequence
    trans_table = int(hit.rsplit('|',1)[-1])
    edited_name = seq_name
    
    # Carry out pairwise alignment of query sequence with its top hit
    alignments = aligner.align(ref_seq, qry_seq) 
    ref_aln, qry_aln = alignments[0]
    
    if not ('-' in ref_aln) and not ('-' in qry_aln):
        # exact match, nothing to correct.
        qry_aln_edited = qry_seq
        qry_aln_edited_nogaps = qry_seq
        no_issues = True
        output_no_change_count += 1
        
    else: 
        # At least one gap is present in the alignment
        ref_aln_endfix, qry_aln_endfix, endfixed, startgapcounter, endgapcounter = end_correction(ref_aln, qry_aln)
        
        if endfixed == True:
            edited_name = seq_name + '|endfixed'
            endfixed_count += 1
           
        qry_aln_edited = qry_aln_endfix

        if ('-' in ref_aln_endfix) or ('-' in qry_aln_endfix):
            # Scan for gaps in reference (extra bases in query)
            no_issues = False
            gap_runs_ref = []
            current_run = ['r']
            
            for i, base in enumerate(ref_aln_endfix):
                if base == "-":
                    current_run.append(i)
                else:
                    if len(current_run)>1:
                        gap_runs_ref.append(current_run)
                        current_run = ['r']
            if len(current_run)>1:
                gap_runs_ref.append(current_run)
               
            # Scan for gaps in query (missing bases)
            gap_runs_qry = []
            current_run = ['q']
            
            for i, base in enumerate(qry_aln_endfix):
                if base == "-":
                    current_run.append(i)
                else:
                    if len(current_run)>1:
                        gap_runs_qry.append(current_run)
                        current_run = ['q']
            if len(current_run)>1:
                gap_runs_qry.append(current_run)

            gap_runs_combined = []

            for q in gap_runs_qry:
                gap_runs_combined.append(q)
            for r in gap_runs_ref:
                gap_runs_combined.append(r)
            
            gap_runs_combined.sort(key=lambda x: x[1], reverse=True) 
            
            for run in gap_runs_combined:
                same_letter = None
                deletion_index_list = []
                Ns_deleted = 0
                query_or_reference = run.pop(0)
                # Take the first letter of the list to tell us if this is a gap in the query or the reference.
                # This also removes the letter from the current gap run list.
                
                if len(run) <= 2: # only try to correct 1 or 2 consecutive gaps
                    N_index_list = []
                    # First, if the gap is in the reference, are the exact position(s) corresponding to the gap(s) Ns? These are safe to delete even without checking for homopolymers.
                    if query_or_reference == 'r':
                        
                        for gap_index in run:
                            if qry_aln_edited[gap_index] == 'N':
                                qry_aln_edited = qry_aln_edited[:gap_index] + '-' + qry_aln_edited[gap_index+1:]
                                Ns_deleted += 1
                                N_index_list.append(gap_index)
                        if Ns_deleted > 0:
                            print('Deleted ambiguous nucleotides from sequence', seq_name,'at position(s)', N_index_list)
                            autocorrect_bool_list.append(True)
    
                        if Ns_deleted == len(run):
                            # We've made all the corrections we need by deleting Ns. Move on to the next gap run or sequence
                            continue
                    
                    # Otherwise, do the homopolymer check.
                    letter, placement, homopolymer, left_start_index, right_end_index = homopolymer_length(qry_aln_edited,run[0])
                    # we are always looking at whether the corresponding position in the QUERY sequences is part of a HP

                    if query_or_reference == 'r':
                        # Gap is in the reference sequence.
                        # Is the corresponding position in the query sequence part of a homopolymer?
                        # If we are deleting more than 1, are both/all letters the same?
                        # If both true, delete that number of letters from the query (replace with gaps)
           
                        same_letter = True
                        if len(homopolymer) >= 4:                           
                            # First, are there any Ns in the homopolymer? If so, preferentially delete those. 

                            if 'N' in homopolymer:                               
                                updated_hp = homopolymer   
                     
                                for index, character in enumerate(homopolymer):
                                                                        
                                    if (character == 'N') & (Ns_deleted < len(run)):
                                        updated_hp = updated_hp[:index] + '-' + updated_hp[index+1:]
                                        N_index_list.append(index+left_start_index)
                                        Ns_deleted += 1
                                
                                # Logic here is different depending on whether current gap is "part" or "adjacent" to the HP
                                
                                if placement == 'part':                                       
                                    
                                    # If deleting on the left
                                    if left_start_index < run[0]:
                                        qry_aln_edited = qry_aln_edited[:left_start_index] + updated_hp + qry_aln_edited[run[-1]:]
                                    
                                    # If deleting on the right
                                    else:
                                        qry_aln_edited = qry_aln_edited[:run[0]] + updated_hp + qry_aln_edited[right_end_index+1:]

                                elif placement == 'adjacent-L':
                                    qry_aln_edited = qry_aln_edited[:left_start_index] + updated_hp + qry_aln_edited[run[0]:]

                                elif placement == 'adjacent-R':
                                    qry_aln_edited = qry_aln_edited[:run[-1]+1] + updated_hp + qry_aln_edited[right_end_index+1:]

                                print('Deleted ambiguous nucleotides from sequence', seq_name,'at position(s)', N_index_list, '(homo)', placement)
                                autocorrect_bool_list.append(True)

                            if Ns_deleted == len(run):
                                # We've made all the corrections we need by deleting Ns. Move on to the next gap run or sequence
                                continue
                            
                            else:                                   
                                if placement == 'part':
                                # The position we are looking at is part of the homopolymer.
                                    for b in run:
                                        if qry_aln_edited[b] != letter:
                                            same_letter = False
                                    if same_letter == True:
                                        qry_aln_edited = qry_aln_edited[:run[0]] + ('-' * len(run)) + qry_aln_edited[(run[-1]+1):]
                                        print('Deleted bases from sequence', seq_name,'at position(s)', run, placement)
                                        autocorrect_bool_list.append(True)
                                    else:
                                        # If the letters in the run are different, we need to instead delete letters from the neighbouring homopolymer.
                                        if left_start_index < run[0]: # homopolymer on the left
                                            qry_aln_edited = qry_aln_edited[:run[0]-len(run)] + ('-' * len(run)) + qry_aln_edited[run[0]:]
                                            deletion_index_list = [x - len(run) for x in run]
                                            
                                        else: # homopolymer on the right
                                            qry_aln_edited = qry_aln_edited[:run[-1]+1] + ('-' * len(run)) + qry_aln_edited[run[-1]+len(run):]
                                            deletion_index_list = [x + len(run) for x in run]

                                        print('Deleted bases from sequence', seq_name,'at position(s)', deletion_index_list, placement)
                                       
                                        autocorrect_bool_list.append(True)
                                        
                                # The homopolymer is adjacent (either to the left or right) of the position we are looking at.
                                # If left, delete the right-most bases. If right, delete the left-most.
                                elif placement == 'adjacent-L':
                                
                                    for d in run:
                                        deletion_index_list.append(d - len(run))
                                            
                                elif placement == 'adjacent-R':
                                    for d in run:
                                         deletion_index_list.append(d + len(run))       
                                            
                                    qry_aln_edited = qry_aln_edited[:deletion_index_list[0]] + ('-' * len(run)) + qry_aln_edited[(deletion_index_list[-1]+1):]
                                    print('Deleted bases from sequence', seq_name,'at position(s)', deletion_index_list, placement)
                                    autocorrect_bool_list.append(True)
                                                                            
                        else:
                            # There's no homopolymer or the length of homopolymer is less than 4.
                            # Only count non-codon indels, i.e. where the length of the gap run is not a multiple of 3.
                            if len(run)%3 != 0: 
                                indels_not_corrected_list.append([query_or_reference, *run])
                                autocorrect_bool_list.append(False)
            
                    elif query_or_reference == 'q':
                        # Gap is in the query sequence.
                        # is this position (in the query) part of a homopolymer?
                        # If so, replace these gaps with Ns, and leave the reference alone.
                        # Query gaps can only be "part" of homopolymer, not "adjacent".
                        if len(homopolymer) >= 4:
                            qry_aln_edited = qry_aln_edited[:run[0]] + ('n' * len(run)) + qry_aln_edited[(run[-1]+1):]
                            print('Added Ns to sequence', seq_name,'at position(s)', run)
                            autocorrect_bool_list.append(True)
                                                        
                        else:
                            # There's no homopolymer or the length of homopolymer is less than 4.
                            # Only count non-codon indels, i.e. where the length of the gap run is not a multiple of 3.
                            if len(run)%3 != 0: 
                                indels_not_corrected_list.append([query_or_reference, *run])
                                autocorrect_bool_list.append(False)

                else:
                    # this corresponds to the length of run if statement. i.e. more than 2 consecutive gaps
                    # Only count non-codon indels, i.e. where the length of the gap run is not a multiple of 3.
                    if len(run)%3 != 0: 
                        indels_not_corrected_list.append([query_or_reference, *run])
                        autocorrect_bool_list.append(False)
                    
        else:
            # this corresponds to the gap if statement. i.e. there are no gaps.
            no_issues = True
            if endfixed == True:
                no_changes_except_endfix_count += 1
            else:    
                output_no_change_count += 1 
    # we are now back in the main loop.

    if True in autocorrect_bool_list:
        # we made at least one autocorrection
        edited_name = edited_name + '|autoedit'
        autocorrect_count += 1
        
    if False in autocorrect_bool_list:
        # Non-codon indels are present which were not autocorrected. 
        # this is where we need the second tier check; trying to correct based only on alignment with the reference.
        
        #First, are there HMM issues? If not, we don't care, make no changes.
        stop_codon, length_issue = hmm_check(qry_aln_edited, trans_table)

        if stop_codon == False and length_issue == False:
        # No HMM issues
            indel_not_corrected_count += 1
            hmm_check_needed = False
            edited_name = edited_name + '|indel'
            
        else:
        # There is a stop codon, and/or the overall length(-1) is not a multiple of 3
            qry_aln_second_tier = second_tier_check(seq_name, qry_aln_edited, indels_not_corrected_list)
            stop_codon_second, length_issue_second = hmm_check(qry_aln_second_tier, trans_table)

            if stop_codon_second == True or length_issue_second == True:
                # There are still stops and/or length issues even after all the edits. Discard the second-tier edits and return the original sequence
                indel_not_corrected_count += 1
                edited_name = edited_name + '|indel'
    
            else:
                 
                non_hp_autocorrect_count += 1
                non_hp_edit = True
                temp_name = edited_name + '|autoedit_nonHP'
                edited_name = edited_name + '|check'
                uncertain_dict[temp_name] = qry_aln_second_tier
            hmm_check_needed = False

    if hmm_check_needed == True:
        stop_codon, length_issue = hmm_check(qry_aln_edited, trans_table)

    if stop_codon == True:
        edited_name = edited_name + '|STOP'
        stop_count += 1

    qry_aln_edited_nogaps = qry_aln_edited.replace('-','')

    # After doing all possible edits, now check the final sequence for length.

    final_seq_length = len(qry_aln_edited_nogaps)
    if not(640 <= final_seq_length <= 670):
       # NUMT; exclude based on length.
       print('Sequence', seq_name,'excluded due to length of', final_seq_length,'bp.')
       problem_dict[edited_name] = qry_aln_edited_nogaps
       excluded_for_length_count += 1
       continue

    if no_issues == False or stop_codon == True or length_issue == True:

        qry_aln_edited = ('-' * startgapcounter) + qry_aln_edited + ('-' * endgapcounter)
        seq_name_parts = seq_name.split('|')
        align_file_name = seq_name_parts[0] + '_' + seq_name_parts[1] + '_align.fasta'
        with chdir(alignment_directory):
           with open(align_file_name, "w") as f:
               f.write(f">{hit}|REF\n{ref_aln}\n")              
               if qry_aln_edited_nogaps == qry_seq and non_hp_edit == False:
                   # If no changes were made to the sequence we only need to output 1 sequence.
                   # However the edited name should be used as it is possible it contains a Stop codon.
                   # Also, it must be the edited sequence so any gaps are retained to align correctly
                   f.write(f">{edited_name}\n{qry_aln_edited}\n")
               else:
                   f.write(f">{seq_name}|orig\n{qry_aln}\n")
                   f.write(f">{edited_name}\n{qry_aln_edited}\n")
                   if non_hp_edit == True:
                   # If there is a sequence with suggested edits based on the second tier check, add this as a final item
                       qry_aln_second_tier = ('-' * startgapcounter) + qry_aln_second_tier + ('-' * endgapcounter)
                       f.write(f">{temp_name}\n{qry_aln_second_tier}\n")

    output_dict[edited_name] = qry_aln_edited_nogaps

output_file = 'edited_seqs.fasta'
problem_file = 'problem_seqs.fasta'
uncertain_file = 'uncertain_edits.fasta'

with open(problem_file, "w") as f:
    for header, qry_seq_edited in problem_dict.items():
        f.write(f">{header}\n{qry_seq_edited}\n")

with open(uncertain_file, "w") as f:
    for header, qry_seq_edited in uncertain_dict.items():
        f.write(f">{header}\n{qry_seq_edited}\n")

with open(output_file, "w") as f:
    for header, qry_seq_edited in output_dict.items():
        f.write(f">{header}\n{qry_seq_edited}\n")

print()
print(len(qry_seq_list), 'sequences analyzed.')
print()
print(excluded_for_length_count, 'were excluded based on length.')
print(excluded_for_nonmatch_count, 'were excluded due to failure to find a match in the reference library.')
print(output_no_change_count,'were output with no changes or issues.')
print(endfixed_count,'were end-corrected.')
print(no_changes_except_endfix_count, 'were output with no changes or issues other than end correction.')
print(autocorrect_count, 'were automatically corrected based on homopolymers (safe).')
print(non_hp_autocorrect_count, 'have corrections suggested based on alignment only (should be checked.)')
print(indel_not_corrected_count, 'contained indels which could not be automatically corrected.')
print(stop_count, 'contained Stop codons.')
print()
print(len(output_dict), 'final sequences output.')