#!/usr/bin/env python3

"""Functions to analyse genomic regions made of pair of genes."""

# Metadata
__author__ = 'Mainguy Jean - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2019 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.1'
__email__ = 'support.bioinfo.genotoul@inra.fr'
__status__ = 'dev'


from collections import defaultdict
import csv
import os
import numpy as np
import logging
import re
import gzip


def get_gene_pair_info(gene_a, gene_b):
    """
    Get genomic information on the gene pair made of gene_a and gene_b.

    When the sequence is circular, the function compute information on the pair gene_a-->gene_b and gene_b-->gene_a

    :Example:
    >>> gene_a = {'start':50, 'end':'150', 'gene_id':'geneA', 'cog_id':'COG1','strand':'+',}
    >>> gene_a.update({'genome':'Genome1','seqid':'Contig1', 'seqid_length':500,  'circular':'False'})
    >>> gene_b = {'start':200, 'end':'260', 'gene_id':'geneB', 'cog_id':'COG1', 'strand':'+'}
    >>> gene_b.update({'genome':'Genome1','seqid':'Contig1', 'seqid_length':500,  'circular':'False'})
    >>> pair_infos = get_gene_pair_info(gene_a, gene_b)
    >>> pair_infos[0]['distance'] == 211
    True
    """
    dist_infos = []
    gene_1, gene_2 = (gene_a, gene_b) if int(gene_a['start']) <= int(gene_b['start']) else (
        gene_b, gene_a)

    if gene_a['seqid'] != gene_b['seqid']:
        # the gene are not from the same dna mol/scaffold/contig
        return []

    if gene_a['circular'] == 'True' and gene_a['gene_id'] != gene_b['gene_id']:
        # gene a and b seq id is circular ..
        # two distances are then computed
        # gene2 to gene1
        seqid_length = gene_a['seqid_length']
        distance_circ = int(seqid_length) - int(gene_2['start']) + 1 + int(gene_1['end'])
        intergenic_distance_cir=int(seqid_length) - int(gene_2['end']) + 1 + int(gene_1['start'])
        cog_pair = f"{gene_2['cog_id']}|{gene_1['cog_id']}"
        strands = f"{gene_2['strand']}|{gene_1['strand']}"
        genes = f"{gene_2['gene_name']}|{gene_1['gene_name']}"
        ids = f"{gene_2['gene_id']}|{gene_1['gene_id']}"
        positions = f"{gene_2['start']}-{gene_2['end']}|{gene_1['start']}-{gene_1['end']}"
        dist_info_circular = {"cog_pair": cog_pair, "genes": genes, 'ids': ids, "strands": strands, "positions": positions,
                              "distance": distance_circ, 'genome': gene_1['genome'], 'seqid': gene_1['seqid'],
                              "intergenic_distance":intergenic_distance_cir,}
        dist_infos.append(dist_info_circular)

    cog_pair = f"{gene_1['cog_id']}|{gene_2['cog_id']}"
    distance = int(gene_2['end']) - int(gene_1['start']) + 1
    intergenic_distance = int(gene_2['start']) - int(gene_1['end']) + 1

    assert distance > 0, f'Distance {distance} <= 0'
    if distance <= 50:
        warning = f'Very small distance of {distance}nt '
        warning += f'in Genome {gene_1["genome"]}. Sequence: {gene_1["seqid"]}. '
        warning += f'Between Gene1 {gene_1["gene_id"]} and Gene2 {gene_2["gene_id"]}'
        logging.warning(warning)

    if gene_2["gene_id"] == gene_1["gene_id"] and (gene_1['start'], gene_1['end']) != (gene_2['start'], gene_2['end']):
        error_message = f'''Gene1 {gene_1["gene_id"]} and Gene2 {gene_2["gene_id"]} have the same id but different coordinates.
            In Genome {gene_1["genome"]}. Sequence: {gene_1["seqid"]}
            May be caused by possible exons not correctly identified durring gff parsing '''
        raise ValueError(error_message)

    strands = f"{gene_1['strand']}|{gene_2['strand']}"
    genes = f"{gene_1['gene_name']}|{gene_2['gene_name']}"
    ids = f"{gene_1['gene_id']}|{gene_2['gene_id']}"
    positions = f"{gene_1['start']}-{gene_1['end']}|{gene_2['start']}-{gene_2['end']}"

    dist_info = {"cog_pair": cog_pair, "genes": genes, "strands": strands, 'ids': ids,
                 "positions": positions, "distance": distance,
                 "intergenic_distance":intergenic_distance,
                 'genome': gene_1['genome'],
                 'seqid': gene_1['seqid']}

    dist_infos.append(dist_info)

    return dist_infos


def get_region_start_and_end(position_string):
    """
    Get region start and end.

    >> get_region_start_and_end("529951-530454|532135-532674")
    (529951, 532674)
    """
    region_start = int(position_string.split('|')[0].split('-')[0])
    region_end = int(position_string.split('|')[1].split('-')[1])

    return region_start, region_end


def get_seq_if_from_pair_info(pair):
    region_start, region_end = get_region_start_and_end(pair['positions'])
    seq_id = f'{pair["genome"]}|{region_start}-{region_end}'
    return seq_id


def remove_circular_duplicate(pairs):
    """
    When the genome is circular, the same two genes are paired twice. The pair is then duplicated.

    For a duplicated pair made of the same two genes, the version of the pair with the smallest distance is kept
    rule to remove circular duplicate: keep the pair that have the smallest distance

    :Example:

    >>> pair1 = {'ids':'gene1|gene2', 'distance':"100", 'seqid':'scaffold10', genome:'genome1'}
    >>> pair2 = {'ids':'gene2|gene1', 'distance':"100000", 'seqid':'scaffold10', genome:'genome1'}
    >>> pair3 = {'ids':'gene3|gene4', 'distance':"1100", 'seqid':'scaffold01', genome:'genome2'}
    >>> remove_circular_duplicate([pair1, pair2]) == [pair1]
    True
    >>> remove_circular_duplicate([pair1, pair2, pair3]) == [pair1, pair3]
    True

    """
    pair_info_without_duplicate_dict = {}
    removed_circ_duplicate = []
    for p in pairs:
        id_pair = tuple([p['seqid'], p['genome']] + sorted(p["ids"].split('|')))
        if id_pair not in pair_info_without_duplicate_dict:
            pair_info_without_duplicate_dict[id_pair] = p
        elif int(p['distance']) < int(pair_info_without_duplicate_dict[id_pair]['distance']):
            removed_circ_duplicate.append(pair_info_without_duplicate_dict[id_pair])
            pair_info_without_duplicate_dict[id_pair] = p
        else:
            removed_circ_duplicate.append(p)
    pair_info_without_duplicate = [p for p in pair_info_without_duplicate_dict.values()]
    logging.info(f'{len(removed_circ_duplicate)} duplicated pairs have been removed')
    return pair_info_without_duplicate


def get_genes_by_genomes(pair_infos):
    """Get genes by genomes."""
    genes_by_genomes = defaultdict(set)
    for p in pair_infos:
        genes_by_genomes[p['genome']] |= set(p['ids'].split('|'))
    return genes_by_genomes


def identify_duplicated_pairs(pair_infos, distance_median):
    """Identify duplicated pairs."""
    len_before_fitering = len(pair_infos)
    pairs_by_genomes = defaultdict(list)
    alternative_pairs = []
    for d in pair_infos:
        pairs_by_genomes[d['genome']].append(d)

    for g, l in pairs_by_genomes.items():
        sorted_pairs_by_diff_median = sorted(
            l, key=lambda k: abs(int(k['distance'])-distance_median))
        for alternative_p in sorted_pairs_by_diff_median[1:]:
            alternative_pairs.append(alternative_p)
            pair_infos.remove(alternative_p)
            # input()

    assert len_before_fitering == len(alternative_pairs) + len(pair_infos)
    return alternative_pairs


def extrem_std_management(distances):
    """
    distances has to be already filtered
    alternative pair have been removed: we have one distance by genome
    Otherwise the final nb associate dwith a std cutoff is the number of pair that give a std < cutoff
    """
    median = np.median(distances)
    distances_sort = sorted(distances, key=lambda x: abs(x-median))
    # print(distances_sort)
    std_variation = [np.std(distances_sort[:i+1]) for i in range(len(distances_sort))]
    # print('std var', std_variation)
    # Get the nb of genome that give a std < the threshold
    std_cutoffs = [100, 500, 1000]
    iter_cutoffs = iter(sorted(std_cutoffs))
    std_cutoff_nb_genomes = {}

    cutoff = next(iter_cutoffs)
    for i, std in enumerate(std_variation):
        while std > cutoff:
            # print(f'std {std} > cutoff {cutoff}: {i}')
            # i start at 0
            # when std > cutoff --> i = number of genome that give std < cutoff
            std_cutoff_nb_genomes[f'nb_genome_std_{cutoff}'] = i

            try:
                cutoff = next(iter_cutoffs)
            except StopIteration:
                break

    # in case of std with all distance < threshold
    # we need to had the key to dict
    # print(std, cutoff)

    if std <= cutoff:
        # print(f'std {std} < cutoff {cutoff}')
        std_cutoff_nb_genomes[f'nb_genome_std_{cutoff}'] = len(std_variation)
        for cutoff in iter_cutoffs:
            # print(f'remaining cutoffs : {cutoff}')
            std_cutoff_nb_genomes[f'nb_genome_std_{cutoff}'] = len(std_variation)

    # print(f'{std_cutoff_nb_genomes}')
    return std_cutoff_nb_genomes, std_variation


def get_pair_organisation(dist_info):
    """
    Get the organisation of the pair

    :Example:

    >>> pair1 = {'cog_pair':'cogA|cogB', 'strands':"+|+"}
    >>> get_pair_organisation(pair1)
    'cogA>-cogB>'

    >>> pair2 = {'cog_pair':'cogB|cogA', 'strands':"-|+"}
    >>> get_pair_organisation(pair2)
    '<cogA-cogB>'

    >>> pair3 = {'cog_pair':'cogB|cogA', 'strands':"-|-"}
    >>> get_pair_organisation(pair3)
    'cogA>-cogB>'
    """
    # if dist_info['distance'] is "":
    #     raise valueError('Distance for the pair is an empty string')
    # cog1 is before cog2 on the plus strand of the sequence
    cog1, cog2 = dist_info['cog_pair'].split('|')

    # sort COG to define cog A and cog B in the pair
    cogA, cogB = sorted((cog1, cog2))

    # strand1 is strand of cog1 and cog2 is strand of cog2
    strands = dist_info['strands']
    strand1, strand2 = dist_info['strands'].split('|')

    if (cog1 == cogA and strands == '+|+') or (cog2 == cogA and strands == '-|-'):
        # cogA FOLOWS cogB
        # A to B --=COG1=(A)=>---=COG2=(B)=>--
        # or
        # -B to -A  --<=COG2=(-B)=---<=COG2=(-A)=--
        organistaion = f'{cogA}>-{cogB}>'

    elif (cog1 == cogB and strands == '+|+') or (cog1 == cogA and strands == '-|-'):
        # cogB FOLOWS cogA
        # B to A --=COG1=(B)=>---=COG2=(A)=>--
        # or
        # -A to -B  --<=COG2=(-B)=---<=COG2=(-B)=--
        organistaion = f'{cogB}>-{cogA}>'

    elif strands == '-|+':
        # --<=cog1=(AorB)=---=cog2=(AorB)=>
        organistaion = f'<{cogA}-{cogB}>'

    elif strands == '+|-':
        # --=cog1=(AorB)=>---<=cog2=(AorB)=
        organistaion = f'{cogA}>-<{cogB}'
    else:
        raise valueError(f'No organisation have been determined for pair: {dist_info}')
    return organistaion


def analyse_cog_pair_old(cog_pair_file):
    """

    """
    organisations = defaultdict(int)
    pair_infos = tools.tsv_to_list_of_dict(cog_pair_file)

    nb_genomes = len({d['genome'] for d in pair_infos})

    pair_infos = remove_circular_duplicate(pair_infos)

    for pair_info in pair_infos:
        orga = get_pair_organisation(pair_info)
        pair_info['organisation'] = orga
        organisations[orga] += 1

    # orga_without_none = {o: c for o, c in organisations.items() if o != None}

    main_orga = max(organisations, key=organisations.get)
    # print(organisations)
    # print(main_orga)
    pairs_main_orga = [l for l in pair_infos if l['organisation'] == main_orga]

    # print([d['distance'] for d in pairs_main_orga])
    median_pre_filtering = np.median([int(d['distance']) for d in pairs_main_orga])
    # print('median pre filtering', median_pre_filtering)

    alternative_pairs = identify_duplicated_pairs(pairs_main_orga, int(median_pre_filtering))

    distances = [int(d['distance']) for d in pairs_main_orga]
    median_post_filtering = np.median(distances)
    # Problematic pair: is a duplicated pair that has a similar distance
    # and orientation than the other pairs
    # All the alternative pair are in the same orientation
    # So only the distance is checked here
    # if the pair distance is close to the median with a factor 10
    # then we consider the pair as problematic
    problematic_pairs = [p for p in alternative_pairs if int(p['distance'])
                         < median_post_filtering*10 and int(p['distance']) > median_post_filtering/10]

    # Alternative genes that would duplicate primer site in the genomes
    genes_by_genomes = get_genes_by_genomes(pair_infos)
    # when the pair is made of the same COG
    # for example COG0093 vs COG0093
    # the computation with the number of alternative genes is affected
    nb_cog_diff_in_pair = len(set(pair_infos[0]['cog_pair'].split('|')))

    nb_alternative_genes = sum(
        [len(genes) - nb_cog_diff_in_pair for genes in genes_by_genomes.values()])

    distances = [int(d['distance']) for d in pairs_main_orga]

    std_thresholds_genome_prct, std_variation = extrem_std_management(distances)

    # and cog_pair[0] != cog_pair[1]:
    # if np.std(distances) < 1000 and len(set(pairs_main_orga[0]['cog_pair'].split('|'))) != 1:
    pair_summary = {'cog_pair': os.path.basename(cog_pair_file[:cog_pair_file.rindex('.')]),
                    "nb_represented_genomes": nb_genomes,
                    "nb_of_pairs": len(pair_infos),
                    "organisations": '|'.join([f'{k}:{v}' for k, v in organisations.items()]),
                    "nb_pair_with_main_orga": len(pairs_main_orga) + len(alternative_pairs),
                    "nb_genome_with_main_orga": len({d['genome'] for d in pairs_main_orga}),
                    'min':  min(distances),
                    'max':  max(distances),
                    'mean': np.mean(distances),
                    'std':  np.std(distances),
                    'median':  np.median(distances),
                    "nb_alternative_pairs": len(alternative_pairs),
                    "nb_problematic_pairs": len(problematic_pairs),
                    "nb_alternative_sites": nb_alternative_genes,
                    'file': cog_pair_file}
    pair_summary.update(std_thresholds_genome_prct)

    return pair_summary, alternative_pairs, pairs_main_orga


def get_stat_of_cog_found_in_genomes(cog_genomic_position_file):

    capture_assemby_name = re.compile("^([A-Z]{3}_\d+.\d+)")

    protein_count = 0
    genome_count = 0
    single_copy_count = 0
    nb_prot_by_genome = defaultdict(int)

    with open(cog_genomic_position_file) as csvfile:
        # header = next(csvfile)[1:]
        # print(header)
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            protein_count += 1
            # remove unecessery last part of the assembly name
            assembly_name = capture_assemby_name.search(row['genome']).group(1)
            nb_prot_by_genome[assembly_name] += 1

    single_copy_count = [
        nb_prot_by_genome for nb_prot_by_genome in nb_prot_by_genome.values()].count(1)

    summary_count = {'protein_count': protein_count,
                     'genome_count': genome_count,
                     'single_copy_count': single_copy_count}
    genome_set = set(nb_prot_by_genome)
    return summary_count, genome_set


def analyse_cog_pair(cog_pair_file, max_length_cutoff):
    """

    """
    organisations = defaultdict(int)
    all_pairs = tools.tsv_to_list_of_dict(cog_pair_file)
    logging.info(f'{len(all_pairs)} pairs in {cog_pair_file}')

    nb_genomes = len({d['genome'] for d in all_pairs})

    all_pairs = remove_circular_duplicate(all_pairs)

    logging.debug(f'Pair found in {nb_genomes} genomes')
    logging.debug(f'{len(all_pairs)} pairs after remove_circular_duplicate')
    assert len(all_pairs) >= nb_genomes

    for pair_info in all_pairs:
        orga = get_pair_organisation(pair_info)
        pair_info['organisation'] = orga
        organisations[orga] += 1

    main_orga = max(organisations, key=organisations.get)

    pairs_main_orga = [p for p in all_pairs if p['organisation'] == main_orga]
    pairs_alternative_orga = [p for p in all_pairs if p['organisation'] != main_orga]

    pairs_main_orga_filtered_length = [
        l for l in pairs_main_orga if int(l['distance']) <= max_length_cutoff]
    pairs_main_orga_exceeding_length = [
        l for l in pairs_main_orga if int(l['distance']) > max_length_cutoff]

    if len(pairs_main_orga_filtered_length) > 0:
        median_pre_filtering = np.median([int(d['distance'])
                                          for d in pairs_main_orga_filtered_length])
        # duplicated pais with main orga and with length < cutoff
        problematic_pairs = identify_duplicated_pairs(
            pairs_main_orga_filtered_length, int(median_pre_filtering))

        distances = [int(d['distance']) for d in pairs_main_orga_filtered_length]

        median = np.median(distances)
        std = np.std(distances)
        mean = np.mean(distances)
        mini = min(distances)
        maxi = max(distances)
    else:
        problematic_pairs = []
        median, std, mean, mini, maxi = None, None, None, None, None
    # Alternative genes that would duplicate primer site in the genomes
    genes_by_genomes = get_genes_by_genomes(all_pairs)

    # when the pair is made of the same COG
    # for example COG0093 vs COG0093
    # the computation with the number of alternative genes is affected
    nb_cog_diff_in_pair = len(set(all_pairs[0]['cog_pair'].split('|')))

    nb_alternative_genes = sum(
        [len(genes) - nb_cog_diff_in_pair for genes in genes_by_genomes.values()])

    # std_thresholds_genome_prct, std_variation = extrem_std_management(distances)

    # and cog_pair[0] != cog_pair[1]:
    # if np.std(distances) < 1000 and len(set(pairs_main_orga[0]['cog_pair'].split('|'))) != 1:
    cog_pair = os.path.basename(cog_pair_file[:cog_pair_file.rindex('.')])

    pair_details = {'cog_pair': cog_pair,
                    'all_pairs': all_pairs,
                    "pairs_main_orga": pairs_main_orga,
                    "pairs_alternative_orga": pairs_alternative_orga,
                    "pairs_main_orga_exceeding_length": pairs_main_orga_exceeding_length,
                    "problematic_pairs_post_len_filtering": problematic_pairs,
                    "pairs_main_orga_post_len_filtering": pairs_main_orga_filtered_length,
                    }

    pair_summary = {'cog_pair': cog_pair,
                    "nb_represented_genomes": nb_genomes,
                    "nb_of_pairs": len(all_pairs),
                    "organisations": '|'.join([f'{k}:{v}' for k, v in organisations.items()]),
                    "main_orga": main_orga,
                    "nb_pair_with_main_orga": len(pairs_main_orga),
                    "nb_pair_with_alternative_orga": len(pairs_alternative_orga),
                    "nb_genome_with_main_orga": len({d['genome'] for d in pairs_main_orga}),
                    "median_distance_pair_with_main_orga": np.median([int(d['distance']) for d in pairs_main_orga_filtered_length]),
                    "nb_genome_with_alternative_orga": len({d['genome'] for d in pairs_alternative_orga}),
                    "max_length_cutoff": max_length_cutoff,
                    "nb_pairs_main_orga_exceeding_length": len(pairs_main_orga_exceeding_length),
                    "nb_genome_with_main_orga_exceeding_length": len({d['genome'] for d in pairs_main_orga_exceeding_length}),
                    "nb_problematic_pairs_post_len_filtering": len(problematic_pairs),
                    "nb_pairs_main_orga_post_len_filtering": len(pairs_main_orga_filtered_length),
                    "nb_genome_with_main_orga_post_len_filtering": len({d['genome'] for d in pairs_main_orga_filtered_length}),
                    'min_post_len_filtering':  mini,
                    'max_post_len_filtering':  maxi,
                    'mean_post_len_filtering': mean,
                    'median_post_len_filtering':  median,
                    'std_post_len_filtering':  std,
                    'nb_alternative_site': nb_alternative_genes,
                    }

    return pair_summary, pair_details


def get_genomic_locations(cds_ids_file, gff_file):
    """
    Get genomic locations.

    First extract CDS ids from the protein_ids_file and
    then retrieve info on CDS in the gff file

    :param protein_ids_file: tsv file containing CDS ids
    :param gff_file:
    :type genes_id_file: str
    :type gff_file: str

    :return genes_info: a dictionary with protein ids as key and another dictionary as value with valuable information
     such as location, seqid and strand and cog_id.
    :rtype: dict
    """
    name_to_cog_dict = parse_ids_file(cds_ids_file)

    gene_infos = retrieve_gene_info(name_to_cog_dict, gff_file)

    # Add cog id into gene_info dict
    for gene_info in gene_infos:
        gene_name = gene_info['gene_name']
        gene_info['cog_id'] = name_to_cog_dict[gene_name]

    return gene_infos


def parse_ids_file(ids_file):
    """
    """
    name_to_cog_dict = {}
    with open(ids_file) as csvfile:
        reader = csv.reader(csvfile,  delimiter='\t')
        for row in reader:
            if row[0].startswith('#'):
                continue
            name_to_cog_dict[row[0]] = row[2]

    return name_to_cog_dict


def retrieve_gene_info(seq_names, gff_file):
    """
    Retieve genes info in the gff using protein ids.

    :param seq_names: protein ids
    :type seq_names: list or dict
    :param gff_file: Gff file of the genome from which the protein ids are from
    :type gff_file: str

    :return gene_info: a dictionary with protein ids as key and another dictionary as value with valuable information
     such as location, seqid and strand.
    :rtype: dict

    # :Example:
    #
    # >>> retrieve_gene_info([''], gff_file)
    # 20090

    .. seealso:: chooseGenes()

    .. note:: see gff fields: https://www.ensembl.org/info/website/upload/gff.html
    """
    gene_infos_by_gene_id = {}
    prot_name_found_count = defaultdict(int)
    circular_regions = {}
    length_regions = {}
    name_pattern = re.compile("Name=([^;]+)")
    id_pattern = re.compile("ID=([^;]+)")
    proper_open = gzip.open if gff_file.endswith('.gz') else open
    with proper_open(gff_file, "rt") as fl:
        for l in fl:
            if l.startswith('#'):
                continue

            seqid, source, type, start, end, score, strand, phase, attributes = l.rstrip().split('\t')

            if type == 'region':
                if 'Is_circular=true' in attributes:
                    circular_regions[seqid] = True
                    length_regions[seqid] = end

                else:
                    circular_regions[seqid] = False
                    length_regions[seqid] = end

            if not type == 'CDS':
                continue
            try:
                prot_name = name_pattern.search(attributes).group(1)
            except AttributeError:
                if 'pseudo=true' not in attributes:
                    logging.info(
                        f'No name found in attributes of the gff line and CDS is not a pseudogene: {attributes}')
                continue
            try:
                gene_id = id_pattern.search(attributes).group(1)
            except AttributeError:
                logging.warning(f'No ID found in attributes of the gff line: {attributes}')
                continue

            if prot_name in seq_names:
                gene_info = {"gene_id": gene_id,
                             "gene_name": prot_name,
                             "seqid": seqid,
                             "start": start,
                             "end": end,
                             "strand": strand,
                             'circular': circular_regions[seqid],
                             'seqid_length': length_regions[seqid]}
                if gene_id in gene_infos_by_gene_id:
                    # COnflict... resolve it.. twice the same ID
                    previous_gene_info = gene_infos_by_gene_id[gene_id]
                    logging.warning(f'Gene id is identical in two line of the gff {os.path.realpath(gff_file)}')
                    logging.warning(f'Gene info {previous_gene_info}')
                    logging.warning(f'Gene info {gene_info}')
                    if previous_gene_info['gene_name'] != gene_info['gene_name']:
                        raise ValueError('gene id is identical but gene name is different..')
                    if previous_gene_info['seqid'] != gene_info['seqid']:
                        raise ValueError('gene id is identical but seqid is different..')
                    if previous_gene_info['strand'] != gene_info['strand']:
                        raise ValueError('gene id is identical but strand is different..')
                    if previous_gene_info['start'] == gene_info['start'] and previous_gene_info['end'] == gene_info['end']:
                        raise ValueError(
                            'Two lines in the gff with identical gene id, seqid, gene name and position...')
                    # Use fct to pair genes...
                    gene_info['genome'] = gff_file
                    previous_gene_info['genome'] = gff_file
                    gene_info['cog_id'] = ''
                    previous_gene_info['cog_id'] = ''
                    previous_gene_info['gene_id'] = 'previous_' + gene_id
                    pair_infos = get_gene_pair_info(gene_info, previous_gene_info)
                    logging.warning(f'Gene with same gene id have been paired: { pair_infos}')
                    # when seqid is circular two possible combination of the two genes exist
                    # Take the combination with the smallest length
                    pair_infos = remove_circular_duplicate(pair_infos)
                    if len(pair_infos) != 1:
                        raise ValueError(
                            f'The two genes have not been correctly joined... {pair_infos}')
                    else:
                        pair_info = pair_infos.pop()

                    new_start, new_end = get_region_start_and_end(pair_info['positions'])
                    new_gene_info = {"gene_id": gene_id,
                                     "gene_name": prot_name,
                                     "seqid": seqid,
                                     "start": new_start,
                                     "end": new_end,
                                     "strand": strand,
                                     'circular': circular_regions[seqid],
                                     'seqid_length': length_regions[seqid]}
                    gene_infos_by_gene_id[gene_id] = new_gene_info
                    logging.warning(f'Fusion of the two genes: {new_gene_info}')
                else:
                    gene_infos_by_gene_id[gene_id] = gene_info
                    prot_name_found_count[prot_name] += 1

    gene_infos = list(gene_infos_by_gene_id.values())
    # Checking...
    if set(prot_name_found_count) != set(seq_names):
        missing = set(seq_names) - set(prot_name_found_count)
        raise ValueError(
            f'Not all genes have been found in the gff ({os.path.realpath(gff_file)}). Missing genes: {missing}')

    for prot_name, count in prot_name_found_count.items():
        if count > 1:
            logging.warning(
                f'Protein {prot_name} has been found {count} times in assembly {os.path.realpath(gff_file)}.')
    return gene_infos


# def resolve_multiple_gene_id_lines(gene_infos):


if __name__ == '__main__':
    import doctest
    doctest.testmod()
