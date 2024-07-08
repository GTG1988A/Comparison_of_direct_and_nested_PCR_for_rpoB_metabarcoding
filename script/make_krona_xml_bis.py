#!/usr/bin/env python3

"""
Make an xml formated to be read by the krona tool ktImportXML

:Example:
$ cat formated_taxo.tsv | python make_krona_xml.py  --color_label 'representative_genome' > test.xml && ktImportXML test.xml -o test.html
$ cat formated_taxo.tsv | python make_krona_xml.py  --color_label 'representative_genome' > test.xml && ktImportXML test.xml -o test.html
$ cat formated_taxo.tsv | python make_krona_xml.py  > test.xml && ktImportXML test.xml -o test.html

"""

# Metadata
__author__ = 'Mainguy Jean - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2019 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.1'
__email__ = 'support.bioinfo.genotoul@inra.fr'
__status__ = 'dev'


import xml.etree.ElementTree as ET
from xml.dom import minidom
import csv
from collections import defaultdict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
import logging
import sys


def build_nested_tree(nested_dict,  taxonomy):

    try:
        taxon = next(taxonomy)
    except StopIteration:
        return

    if taxon in nested_dict:
        next_nested_dict = nested_dict[taxon]
    else:
        next_nested_dict = {}
        nested_dict[taxon] = next_nested_dict

    build_nested_tree(next_nested_dict, taxonomy)


def process_labeled_taxonomies(tsvfile):
    nested_taxo_dict = {}

    reader = csv.DictReader(tsvfile, delimiter='\t')

    taxon_count_by_label = {label.replace(' ', '_'): defaultdict(int)
                            for label in reader.fieldnames if label != "taxonomy"}

    taxon_count_by_label['database'] = defaultdict(int)
    for line in reader:
        taxonomy_line = line["taxonomy"].rstrip().split(';')[::-1]
        labels_dict = [k for k, v in line.items() if k !=
                       'taxonomy' and (v == '1' or v == 'True')]

        for t in taxonomy_line:
            taxon_count_by_label['database'][t] += 1
            for label in labels_dict:
                taxon_count_by_label[label.replace(' ', '_')][t] += 1
        build_nested_tree(nested_taxo_dict, iter(taxonomy_line))

    return nested_taxo_dict, taxon_count_by_label


def process_labeled_taxonomies_test(tsvfile, dataset_labels, color_label=None):
    nested_taxo_dict = {}

    reader = csv.DictReader(tsvfile, delimiter='\t')
    # taxon_count_by_label = {label.replace(' ', '_'): defaultdict(int)
    #                         for label in reader.fieldnames if label != "taxonomy" and not label.stratswith(color_label)}
    taxon_count_by_label = {}
    # dataset_labels = []
    taxon_count_by_label = {label.replace(' ', '_'): defaultdict(int)
                            for label in reader.fieldnames if label != "taxonomy" and label not in dataset_labels and label != color_label}
    logging.info(f'in process {taxon_count_by_label}')
    if dataset_labels != [l for l in dataset_labels if l in reader.fieldnames]:
        raise NameError('dataset label given does not correspond to col name of the tax file')

    taxon_count_by_label['database'] = {}
    if color_label:
        taxon_count_by_label[color_label] = {}
    for line in reader:
        taxonomy_line = line["taxonomy"].rstrip().split(';')[::-1]
        col_dict = [k for k, v in line.items() if k != 'taxonomy' and (v == '1' or v == 'True')]

        for t in taxonomy_line:
            # taxon_count_by_label['database'][t] += 1
            if t not in taxon_count_by_label['database']:
                taxon_count_by_label['database'][t] = {col: 0 for col in dataset_labels}
                if color_label:
                    taxon_count_by_label[color_label][t] = {col: 0 for col in dataset_labels}
            for col_name in col_dict:
                if col_name in dataset_labels:
                    taxon_count_by_label['database'][t][col_name] += 1
                    if color_label in col_dict:
                        taxon_count_by_label[color_label][t][col_name] += 1
                elif col_name != color_label:
                    taxon_count_by_label[col_name.replace(' ', '_')][t] += 1

        build_nested_tree(nested_taxo_dict, iter(taxonomy_line))
    logging.info(f'last {taxon_count_by_label.keys()}')
    return nested_taxo_dict, taxon_count_by_label


def pretty(d, indent=0):
    for key, value in d.items():
        print('\t' * indent + str(key))
        if isinstance(value, dict):
            pretty(value, indent+1)
        else:
            print('\t' * (indent+1) + str(value))


def manage_nodes(childrens, parent, taxon_count_by_label, nb_datasets=1):
    for child in childrens:
        # add child to xml obj
        # print(type(parent))
        child_element = ET.SubElement(parent, 'node', attrib={'name': child})
        for label in taxon_count_by_label:
            label_element = ET.SubElement(child_element, label)
            if type(taxon_count_by_label[label][child]) == dict:
                for dataset_label, count in taxon_count_by_label[label][child].items():
                    val = ET.SubElement(label_element, 'val')
                    val.text = str(count)
            else:
                for _ in range(nb_datasets):
                    val = ET.SubElement(label_element, 'val')
                    val.text = str(taxon_count_by_label[label][child])

        next_childrens = childrens[child]
        manage_nodes(next_childrens, child_element, taxon_count_by_label, nb_datasets)


def parse_arguments():
    parser = ArgumentParser(description="...",
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("--name", type=str, default='',
                        help="prefixed to be displayed before the label (can be the name of the region ie rpoB)")
    parser.add_argument("--color_label", type=str, default=None,
                        help="Label to which apply coloring")
    parser.add_argument("--dataset_labels", nargs='+', type=str, default=[],
                        help="Label to which apply coloring")
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                        action="store_true")
    parser.add_argument('taxonomy_labeled_file', nargs='?', type=FileType('r'),
                        default=sys.stdin)
    parser.add_argument('xml_outfile', nargs='?', type=FileType('w'),
                        default=sys.stdout)
    parser.add_argument('--color_by_default', help="color krona by default",
                        action="store_true", default=False)

    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info('Mode verbose ON')

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    taxonomy_fl = args.taxonomy_labeled_file
    region_name = args.name
    xml_outfile = args.xml_outfile
    color_label = args.color_label
    dataset_labels = args.dataset_labels
    logging.info(dataset_labels)
    logging.info(taxonomy_fl)
    # nested_taxo_dict, taxon_count_by_label = process_labeled_taxonomies(taxonomy_fl)

    nested_taxo_dict, taxon_count_by_label = process_labeled_taxonomies_test(
        taxonomy_fl, dataset_labels, color_label)
    logging.info(f'dataset_labels {len(dataset_labels)} taxon_count_by_label {taxon_count_by_label.keys()}')
    # print(taxon_count_by_label.keys())
    krona = ET.Element('krona')

    attributes = ET.SubElement(krona, 'attributes', attrib={'magnitude': f'database'})

    attribute = ET.SubElement(attributes, 'attribute', attrib={'display': "Number of genomes"})
    attribute.text = 'database'
    labels = [label for label in taxon_count_by_label if label !=
              'database' and label not in dataset_labels]
    logging.info(f'LABELS not in dataset labels {labels}')
    for label in labels:  # if more than just the database info...
        attribute = ET.SubElement(attributes, 'attribute', attrib={'display': f'{region_name} percentage {label}'})
        attribute.text = f'{label}_prct'

        attribute = ET.SubElement(attributes, 'attribute', attrib={'display': f'{region_name} {label}'})
        attribute.text = f'{label}'
        if label == color_label:
            taxon_count_by_label[f'{color_label}_prct'] = {}
            for t, count_in_datasets in taxon_count_by_label['database'].items():
                if t not in taxon_count_by_label[f'{color_label}_prct']:
                    taxon_count_by_label[f'{color_label}_prct'][t] = {col: 0 for col in dataset_labels}
                for dataset_label, total_count in count_in_datasets.items():
                    count_color = taxon_count_by_label[color_label][t][dataset_label]
                    if total_count == 0:
                        taxon_count_by_label[f'{color_label}_prct'][t][dataset_label] = 0
                    else:
                        taxon_count_by_label[f'{color_label}_prct'][t][dataset_label] = round(100*(count_color/total_count), 2)
        else:
            taxon_count_by_label[f'{label}_prct'] = {t: round(100 * taxon_count_by_label[label][t] /
                                                              taxon_count_by_label['database'][t], 2) for t in taxon_count_by_label['database']}

    logging.info(taxon_count_by_label.keys())

    datasets = ET.SubElement(krona, 'datasets')
    for dataset_label in dataset_labels:
        dataset = ET.SubElement(datasets, 'dataset')
        dataset.text = dataset_label

    if color_label:
        if color_label in labels:
            default_color = 'true' if args.color_by_default else 'false'
            color = ET.SubElement(krona, 'color', attrib={"attribute": f'{color_label}_prct',
                                                          "valueStart": "0", "valueEnd": "100",
                                                          "hueStart": "0", "hueEnd": "120", 'default': default_color})
            color.text = ' '
        else:
            logging.warning(f'Provided color label "{color_label}" is not found in the columns of the input file that can be colored : {labels}')

    # f'{label}_coverage': taxon_count_by_label[f'{label}_coverage']}
    logging.info(f'dataset_labels {len(dataset_labels)} taxon_count_by_label {taxon_count_by_label.keys()}')
    # print(minidom.parseString(ET.tostring(krona)).toprettyxml(indent="   "))
    manage_nodes(nested_taxo_dict, krona, taxon_count_by_label, nb_datasets=len(dataset_labels))
    # pretty(nested_taxo_dict)
    # pretty(taxon_count_by_label)
    mydata = minidom.parseString(ET.tostring(krona)).toprettyxml(indent="   ")

    xml_outfile.write(mydata)
    # with open("items2.xml", "w") as fl:
    #     fl.write(mydata)


if __name__ == '__main__':
    main()
