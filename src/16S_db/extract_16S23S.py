#!/usr/bin/env python3

"""
Retrieve genomic positions and extract 16S sequences

:Example:
python
"""

# Metadata
__author__ = 'Mainguy Jean - Plateforme bioinformatique Toulouse & Agoutin Gabryelle INRAE Auzeville'
__copyright__ = 'Copyright (C) 2022 INRAE'
__license__ = 'GNU General Public License'
__version__ = '0.1'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'dev'


from collections import defaultdict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os
import logging
import re
import identify_rrna_gene
from tqdm import tqdm
import sys
from multiprocessing.pool import Pool
# from concurrent.futures import ThreadPoolExecutor as Pool

# import make_cog_pairs

import regions_analysis_fct
from Bio import SeqIO
import gzip
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv
import subprocess


def load_taxonomic_info(taxonomy_file):
    taxonomic_info = {}
    with open(taxonomy_file, 'r') as file:
        reader = csv.reader(file, delimiter=' ')
        for row in reader:
            assembly_name = row[1]
            taxid = row[2]
            taxonomic_info[assembly_name] = taxid
    return taxonomic_info


def get_gff_and_fna_files(assembly_dir):

    file_found_count = 0
    gff_file = None
    fna_file = None
    print(os.listdir(assembly_dir))
    for data_file in os.listdir(assembly_dir):

        if data_file.endswith('.gff') or data_file.endswith('.gff.gz'):
            gff_file = os.path.join(assembly_dir, data_file)
            file_found_count += 1

        if data_file.endswith('.fna') or data_file.endswith('.fna.gz'):
            fna_file = os.path.join(assembly_dir, data_file)
            file_found_count += 1

    return gff_file, fna_file


def format_gene_for_extraction(gene_info, store_dir):
    """
    Format gene info to be extracted from the genome by the function extract_sequences
    """
    start, end = int(gene_info['start']),  int(gene_info['end'])
    gene_region = {'region_name': gene_info['cog_id'],
                   'seqid': gene_info['seqid'],
                   'start': start,
                   'end': end,
                   'strand': int(f"{gene_info['strand']}1"),
                   'header': f'{gene_info["genome"]}|{gene_info["seqid"]}|{start}-{end}',
                   'description': f"{gene_info['gene_name']} {gene_info['gene_id']}",
                   'store_dir': store_dir
                   }
    return gene_region



def get_seq_of_gene(record, start, end, strand):
    start= int(start)
    end=int(end)

    if strand == '+':
        strand = +1
    elif strand == '-':
        strand = -1
    else:
        strand = None


    if start < end:
        # everything is normal
        location = FeatureLocation(start-1, end)  # minus one to be in zero base

    elif start >= end:
        # Then the region is overlapping the 0 of the chromosome
        location = FeatureLocation(start-1, len(record)) + FeatureLocation(0, end)
        # print("MULPTIPLE LOCATION ", location)
    else:
        logging.critical('Problem in location of seq')
        exit(1)

    feat = SeqFeature(location,  strand=strand)
    return feat.extract(record)



def extract_sequences(fna_file, regions, taxid):

    nb_region_extracted = 0
    proper_open = gzip.open if fna_file.endswith('.gz') else open
    regions_by_seqid = defaultdict(list)

    sequences_str = ''

    for r in regions:
        regions_by_seqid[r['seqid']].append(r)

    try:
        with proper_open(fna_file, 'rt') as handle:

            for record in SeqIO.parse(handle, "fasta"):
                for region in regions_by_seqid[record.id]:
                    try:
                        print(region['start'], region['end'], region['strand'])
                        seq = get_seq_of_gene(record, region['start'], region['end'], region['strand'])
                        print(seq)
                        #seq.id = region['header']
                        #seq.description = region['description']

                        sequences_str += f">{region['seqid']}| taxid={taxid};\n"
                        sequences_str += f"{str(seq.seq)}\n"

                        nb_region_extracted += 1
                    except Exception as e:
                        logging.error("An error occured while opening {fna_file}: {e}")
                        continue
    except Exception as e:
        logging.error(f"An error occured while opening {fna_file}: {e}")
    assert nb_region_extracted == len(
        regions), f'{nb_region_extracted} regions extracted but {len(regions)} expected'

    return sequences_str

def identify_region_in_gff(gff_file, assembly_name,taxid):


    regions_to_extract = []

    gene_infos, gff_parser_issue_recorder = identify_rrna_gene.retrieve_rrna_genes_from_gff(gff_file, assembly_name,taxid)
    return gene_infos, gff_parser_issue_recorder


def parse_arguments():
    parser = ArgumentParser(description="...",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-a', "--assembly_dirs", type=str, help="Directories of assemblies with at least the gff file and fna file. A path to a filename containing one dir per line.")
    parser.add_argument('-t', "--taxonomy_file", type=str, help="File containing taxonomic information.")
    parser.add_argument('-c','--cpus', type=int, default=1,
                        help="Number of thread to use to extract sequences.")
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")

    args = parser.parse_args()
    return args

def identify_and_extract_16S23S_sequences(assembly_dir,taxonomy_file):

    logging.debug(f'Processing {assembly_dir}')
    if not os.path.isdir(assembly_dir.rstrip()):
        logging.warning(f'Given directory does not exist: {assembly_dir}')
        log = {'Assembly dir does not exist':True}

        return None, None, log


    capture_assemby_acc = re.compile("([A-Z]{3}_\d+.\d+)")
    print(capture_assemby_acc)
    parser_log = {"missing_gff_file_in_assembly": False,
                        "missing_fna_file_in_assembly":False}

    gff_file, fna_file = get_gff_and_fna_files(assembly_dir)
    # logging.debug(f'Processing {gff_file} {fna_file}')
    if gff_file is None:
        logging.critical(f"GFF file is missing in the current assembly dir: {assembly_dir}. skipping this assembly")
        parser_log['missing_gff_file_in_assembly'] =True


    if fna_file is None:
        logging.critical(f"Fna file is missing in assembly dir: {assembly_dir}. skipping this assembly")
        parser_log['missing_fna_file_in_assembly'] = True

    if not fna_file or not gff_file:
        return None, None, parser_log

    assembly_accession = capture_assemby_acc.search(os.path.basename(fna_file)).group(1)
    taxonomic_info = load_taxonomic_info(taxonomy_file)
    taxid = taxonomic_info.get(assembly_accession, "Unknown")

    region_infos, gff_parser_log = identify_region_in_gff(gff_file, assembly_accession, taxid)
    try:
        seq_str = extract_sequences(fna_file, region_infos, taxid)
    except Exception as e:
        logging.error(f"An error occured while opening {fna_file}: {e}")
        seq_str=None

    return seq_str, region_infos, gff_parser_log

def main():
    args = parse_arguments()

    level = logging.DEBUG if args.verbose else logging.INFO

    logging.basicConfig(level=level,
                        format='%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info("Command: " + " ".join([arg for arg in sys.argv]))

    if args.verbose:
        logging.debug('Mode verbose ON')

    assembly_dir_file = args.assembly_dirs
    taxonomy_file = args.taxonomy_file
    cpus = args.cpus



    extracted_seq_file = '16S.fna.gz'
    region_infos_file = '16S.tsv.gz'


    rrna_genes_to_fetch = ["16S"]

    logging.info(f'ribosomal RNA genes : {rrna_genes_to_fetch} are being retrieved in genomes')

    region_found_count = 0
    assemblies_with_region_count = 0
    assemblies_count = 0

    gff_parser_log_counter = defaultdict(int)


    logging.info(f'Sequences are written in {extracted_seq_file}')



    with open(assembly_dir_file, 'r') as assembly_dir_fh:
        # generate a generator with assembly dir
        assemblies_dirs = [assembly_dir.rstrip() for assembly_dir in assembly_dir_fh]

    logging.info(f'Extracting {" ".join(rrna_genes_to_fetch)} in {len(assemblies_dirs)} assemblies')


    with gzip.open(extracted_seq_file, 'wt') as out_seq_fh, gzip.open(region_infos_file, 'wt') as out_info_fh:

        logging.info(f'Information on the extracted region are written in {region_infos_file}')
        headers = ['gene_id', 'seqid','product','gene_name', 'start', 'end', 'strand', 'circular', 'seqid_length','partial','genome_name','species_taxid']
        info_writer = csv.DictWriter(out_info_fh, fieldnames=headers, delimiter='\t')
        info_writer.writeheader()


 # Process each assembly directory
        progress_bar = tqdm(total=len(assemblies_dirs), desc="Progress..", unit="assembly")
        for assembly_dir in assemblies_dirs:
            sequences_str, region_infos, gff_parser_log = identify_and_extract_16S23S_sequences(assembly_dir, taxonomy_file)
            assemblies_count += 1

            if sequences_str is not None:
                out_seq_fh.write(sequences_str)
                for region_info in region_infos:
                    info_writer.writerow(region_info)
                region_found_count += len(region_infos)
                assemblies_with_region_count += 1 if region_infos else 0

            for issue, flag in gff_parser_log.items():
                if flag:
                    gff_parser_log_counter[issue] += 1

            # Mettre à jour la barre de progression
            progress_bar.update(1)

        progress_bar.close()



        logging.info(f'{assemblies_count} assemblies have been scanned for 16S sequences')

        logging.info(f'{region_found_count} 16S have been found in {assemblies_with_region_count} assemblies')


        if gff_parser_log_counter:
            logging.info(f'Summing up gff parsing error:')
            for issue, counter in gff_parser_log_counter.items():
                logging.info(f'{issue}: {counter}/{assemblies_count} assemblies')


if __name__ == '__main__':
    main()
