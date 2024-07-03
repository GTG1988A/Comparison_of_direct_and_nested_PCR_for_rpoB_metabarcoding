#!/usr/bin/env python3

"""
Description

:Example:
identify_rrna_gene.py --gff_files GCF_*/*gff* -v
"""

# Metadata
__author__ = 'Mainguy Jean - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2019 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.1'
__email__ = 'support.bioinfo.genotoul@inra.fr'
__status__ = 'dev'


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
import gzip
import re
import os
import csv

def retrieve_rrna_genes_from_gff(gff_file, assembly_name, taxid):
    """

    .. note:: see gff fields: https://www.ensembl.org/info/website/upload/gff.html
    """
    gene_infos = []
    circular_regions = {}
    length_regions = {}

    partial_genes_count = 0
    complete_genes_count = 0

    name_pattern = re.compile('(\d{1,2}S)\w{0,1} (ribosomal |ribosormal ||r)RNA', re.IGNORECASE) # ribosormal is found in this genbank assembly  GCA_000828835.1. Likely a funny silly human mistake
    name_pattern2 = re.compile('ribosomal RNA[-\ ]{1}(\d{1,2}S)', re.IGNORECASE)
    name_pattern3 = re.compile('(\d{1,2}S) (small) subunit ribosomal RNA', re.IGNORECASE)

    name_pattern4 = re.compile('(small) subunit ribosomal RNA', re.IGNORECASE)
    name_pattern5 = re.compile('(SSU) ribosomal RNA', re.IGNORECASE)



    proper_open = gzip.open if gff_file.endswith('.gz') else open

    issue_recorder = {"product_not_found_in_attributes":False,
                      "ID_not_found_in_attributes":False,
                      "rRNA_product_regex_fail":False}
    with proper_open(gff_file, "rt") as fl:
        for l in fl:
            if l.startswith('#'):
                continue
            if l.startswith('>'):
                # gff3 format contain fasta sequences at the end of the file
                # when sequences has been reach we stop
                break

            seqid, source, type, start, end, score, strand, phase, attributes = l.rstrip().split('\t')

            if type == 'region':
                length_regions[seqid] = end

                if 'Is_circular=true' in attributes:
                    circular_regions[seqid] = True
                else:
                    circular_regions[seqid] = False

            if type == 'rRNA':
                attributes_dict = {attr.split('=')[0]: attr.split('=')[1]
                                   for attr in attributes.split(';')}
                try:
                    product = attributes_dict['product']
                except KeyError:
                    logging.warning(
                        f'product is not found in attributes section of the rRNA gff line: "{attributes}" in {gff_file}')
                    issue_recorder["product not found in attributes"] = True
                    if 'gene' in attributes_dict:
                        logging.info(f'the gene attribute is used to replace missing product attribute in {gff_file}')
                        product = attributes_dict['gene']
                        issue_recorder["use gene attributes to replace missing product attribute"] = True

                    else:
                        continue
                try:
                    gene_id = attributes_dict['ID']
                except KeyError:

                    gene_id = f"{seqid}_{start}_{end}"
                    logging.warning(
                        f'ID is not found in attributes section of the rRNA gff line - line: {l}. new id built up:{gene_id}')
                    issue_recorder["ID_not_found_in_attributes"] = True

                if name_pattern.match(product):
                    # logging.info(name_pattern.search(product).group(1))
                    rna_name = name_pattern.match(product).group(1)

                elif name_pattern2.match(product):
                    rna_name = name_pattern2.match(product).group(1)

                elif name_pattern3.match(product):
                    rna_name = name_pattern3.match(product).group(1)

                elif name_pattern4.match(product):
                    rna_name = name_pattern4.match(product).group(1)
                    rna_name = "16S" if rna_name == 'small' else "23S"

                elif name_pattern5.match(product):
                    rna_name = name_pattern5.match(product).group(1)
                    rna_name = "16S" if rna_name == 'SSU' else "23S"

                else:
                    logging.warning(
                        f'regex patterns did not match rRNA name in product section: {product} in {gff_file}')
                    issue_recorder["rRNA product regex fail"] = True
                    continue

                # print(l)  # print to keep all rrna line of interest..
                partial = 'partial=true' in attributes
                if rna_name == "16S":
                    gene_infos.append({"gene_id": gene_id,
                                       "seqid": seqid,
                                       "product": product,
                                       "gene_name": rna_name,
                                       "start": start,
                                       "end": end,
                                       "strand": strand,
                                       'circular': circular_regions[seqid],
                                       'seqid_length': length_regions[seqid],
                                       'partial': partial,
                                       'genome_name': assembly_name,
                                       'species_taxid': taxid})

        logging.debug(
            f'rRNA identification: {os.path.realpath(gff_file)}: {partial_genes_count + complete_genes_count} rrna genes: {complete_genes_count} complete and {partial_genes_count} partial.')
        logging.debug(f"rRNA identification: rRNA gene found: { {info['gene_name'] for info in gene_infos} }")
        return gene_infos, issue_recorder

def write_dict_to_tsv(dict_to_write, handlers_dict, filename):
    """Write dict object to tsv."""
    if filename in handlers_dict:
        handlers_dict[filename].writerow(dict_to_write)
    else:
        fh = open(filename, 'w', newline='')
        writer = csv.DictWriter(fh, fieldnames=list(dict_to_write.keys()), delimiter='\t')
        writer.writeheader()
        writer.writerow(dict_to_write)
        handlers_dict[filename] = writer
        try:
            handlers_dict['open_files'].append(fh)
        except KeyError:
            handlers_dict['open_files'] = []
            handlers_dict['open_files'].append(fh)
