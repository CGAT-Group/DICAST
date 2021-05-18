#!/usr/bin/env python
# coding: utf-8

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

gtf_file = open(input_file)
irfinder_gtf_file = open(output_file, 'w')
new_description = {}
for line in gtf_file:
    if not line.startswith('#'):
        strings = line.split("\t")
    
        if strings[2] == "exon":
            description = strings[-1].split()
            new_description['gene_id'] = description[1]
            new_description['transcript_id'] = description[3]
            new_description['template'] = description[5]
            new_description['gene_exon_number'] = description[7]
            new_description['tr_start'] = description[9]
            new_description['tr_end'] = description[11].strip() + ";"
        
            new_description['Parent'] = description[3]
            new_description['gene_name'] = description[1]
            new_description['gene_biotype'] = '"processed_transcript";'
            new_description['transcript_biotype'] = '"processed_transcript";'
        
            new_description_list = []
            for k,v in new_description.items():
                new_description_list.append(k)
                new_description_list.append(v)
        
            new_description_string = ' '.join(new_description_list)
        
            new_strings = strings[:-1]
            new_strings.append(new_description_string)
            result = '\t'.join(new_strings)
            irfinder_gtf_file.write(result + "\n")
        else:
            irfinder_gtf_file.write(line)
gtf_file.close()
irfinder_gtf_file.close()
