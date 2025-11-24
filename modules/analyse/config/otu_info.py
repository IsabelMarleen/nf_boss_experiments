

# this translates the names as they appear in the downloaded fasta to cleaner names
# and it also associates contigs from the same organism to the same name
otus_clean_names = {
    # 'chr21' : 'hschr21',
    '21' : 'hschr21' #,
    # 'ENA|CP068264|CP068264.2' : 'hschr14'
}

# this is a dict that just contains the chromosomes, no plasmids
# used for mapping against
otus_no_plasmids = {k: v for k, v in otus_clean_names.items() if 'plasmid' not in k}


# genome sizes used to calculate mean coverage
genome_sizes = {
    'hschr21': 45090682 #,
    #'hschr14': 101161492,
}

 
