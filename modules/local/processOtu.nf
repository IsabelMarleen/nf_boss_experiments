process processOtu{
    input:
        val otu_clean_names
        val otu_genome_sizes
    output:
        path 'otu_info.py'
    script:
    fixed_otu_genome_sizes = otu_genome_sizes.replaceAll('"', "'")
    fixed_otu_clean_names = otu_clean_names.replaceAll('"', "'")
        """
        echo "otus_clean_names = $fixed_otu_clean_names
otus_no_plasmids = {k: v for k, v in otus_clean_names.items() if 'plasmid' not in k}
genome_sizes = $fixed_otu_genome_sizes" > otu_info.py
        """

}