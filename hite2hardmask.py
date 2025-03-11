#!/usr/bin/env python
# coding: utf-8

# In[27]:


# import os
# os.chdir(r'C:\Users\WeiChuanzheng\Desktop\hite')
import argparse

# In[28]:


def read_repeat_positions(file):
    repeat_positions = []    
    with open(file, 'r') as f:
        for _ in range(3):
            next(f)
        
        for line in f:
            columns = line.split()
            chrom = columns[4]
            start = min(int(columns[5]), int(columns[6]))
            end = max(int(columns[5]), int(columns[6]))
            repeat_positions.append((chrom, start, end))
            
    repeat_positions.sort(key=lambda x: (x[0], x[1]))

    merged_positions = []
    for chrom, start, end in repeat_positions:
        if not merged_positions or merged_positions[-1][0] != chrom or merged_positions[-1][2] < start - 1:
            merged_positions.append([chrom, start, end])
        else:
            merged_positions[-1][2] = max(merged_positions[-1][2], end)

    merged_positions = [(chrom, start, end) for chrom, start, end in merged_positions]

    return merged_positions


# In[29]:


def mask_fasta_with_repeats(fasta_file, repeat_positions, output_file):
    chrom_seq = {}
    
    with open(fasta_file, 'r') as f_in:
        current_chrom = None
        current_seq = []
        
        for line in f_in:
            if line.startswith(">"):
                if current_chrom:
                    chrom_seq[current_chrom] = ''.join(current_seq)
                current_chrom = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        if current_chrom:
            chrom_seq[current_chrom] = ''.join(current_seq)
    
    repeat_positions_by_chrom = {}
    for chrom, start, end in repeat_positions:
        if chrom not in repeat_positions_by_chrom:
            repeat_positions_by_chrom[chrom] = []
        repeat_positions_by_chrom[chrom].append((start, end))
    
    for chrom, seq in chrom_seq.items():
        if chrom in repeat_positions_by_chrom:
            repeat_positions = repeat_positions_by_chrom[chrom]
            
            regions = []
            prev_end = 0
    
            for start, end in repeat_positions:
                if start - 1 > prev_end:
                    regions.append(('U', prev_end+1, start-1))
                regions.append(('N', start, end))
                prev_end = end
            if prev_end < len(seq):
                regions.append(('U', prev_end + 1, len(seq)))

            final_seq = []
            for region_type, start, end in regions:
                if region_type == 'N':
                    length = end - start + 1
                    final_seq.append("N" * length)
                else:
                    segment = seq[start - 1:end]
                    final_seq.append(segment.upper())

            chrom_seq[chrom] = ''.join(final_seq)

    with open(output_file, 'w') as f_out:
        for chrom, seq in chrom_seq.items():
            f_out.write(f">{chrom}\n")
            for i in range(0, len(seq), 100):
                f_out.write(seq[i:i + 100] + "\n")

def main():
    parser = argparse.ArgumentParser(description="Mask FASTA genome based on repeat positions with hard masking.")
    
    parser.add_argument('-r', '--repeat', type=str, required=True, help="Path to the repeat file (e.g., tmp.out).")
    parser.add_argument('-g', '--genome', type=str, required=True, help="Path to the genome FASTA file (e.g., Chr01.fasta).")
    parser.add_argument('-o', '--out', type=str, required=True, help="Path to output the masked FASTA file (e.g., masked_output.fasta).")
    
    args = parser.parse_args()
    repeat_positions = read_repeat_positions(args.repeat)
    mask_fasta_with_repeats(args.genome, repeat_positions, args.out)

if __name__ == '__main__':
    main()