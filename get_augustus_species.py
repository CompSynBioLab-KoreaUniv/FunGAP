#!/usr/bin/env python3

'''
Get Augustus species argument by taxon name
Last updated: Aug 12, 2020
'''

import re
import sys
from argparse import ArgumentParser
from collections import defaultdict

from Bio import Entrez

# Parameter
SPECIES_DATA = [
    ('aspergillus_fumigatus', 'Aspergillus', 'Eurotiales', 'Eurotiomycetes',
     'Pezizomycotina', 'Ascomycota'),
    ('aspergillus_nidulans', 'Aspergillus', 'Eurotiales', 'Eurotiomycetes',
     'Pezizomycotina', 'Ascomycota'),
    ('aspergillus_oryzae', 'Aspergillus', 'Eurotiales', 'Eurotiomycetes',
     'Pezizomycotina', 'Ascomycota'),
    ('aspergillus_terreus', 'Aspergillus', 'Eurotiales', 'Eurotiomycetes',
     'Pezizomycotina', 'Ascomycota'),
    ('botrytis_cinerea', 'Botrytis', 'Helotiales', 'Leotiomycetes',
     'Pezizomycotina', 'Ascomycota'),
    ('candida_albicans', 'Candida', 'Saccharomycetales', 'Saccharomycetes',
     'Saccharomycotina', 'Ascomycota'),
    ('candida_guilliermondii', 'Candida', 'Saccharomycetales',
     'Saccharomycetes', 'Saccharomycotina', 'Ascomycota'),
    ('candida_tropicalis', 'Candida', 'Saccharomycetales', 'Saccharomycetes',
     'Saccharomycotina', 'Ascomycota'),
    ('chaetomium_globosum', 'Chaetomium', 'Sordariales', 'Sordariomycetes',
     'Pezizomycotina', 'Ascomycota'),
    ('coccidioides_immitis', 'Coccidioides', 'Onygenales', 'Eurotiomycetes',
     'Pezizomycotina', 'Ascomycota'),
    ('coprinus_cinereus', 'Coprinus', 'Agaricales', 'Agaricomycetes',
     'Agaricomycotina', 'Basidiomycota'),
    ('cryptococcus_neoformans_gattii', 'Cryptococcus', 'Tremellales',
     'Tremellomycetes', 'Agaricomycotina', 'Basidiomycota'),
    ('cryptococcus_neoformans_neoformans_B', 'Cryptococcus', 'Tremellales',
     'Tremellomycetes', 'Agaricomycotina', 'Basidiomycota'),
    ('cryptococcus_neoformans_neoformans_JEC21', 'Cryptococcus', 'Tremellales',
     'Tremellomycetes', 'Agaricomycotina', 'Basidiomycota'),
    ('debaryomyces_hansenii', 'Debaryomyces', 'Saccharomycetales',
     'Saccharomycetes', 'Saccharomycotina', 'Ascomycota'),
    ('encephalitozoon_cuniculi_GB', 'Encephalitozoon', '-', '-', '-',
     'Microsporidia'),
    ('eremothecium_gossypii', 'Eremothecium', 'Saccharomycetales',
     'Saccharomycetes', 'Saccharomycotina', 'Ascomycota'),
    ('fusarium_graminearum', 'Fusarium', 'Hypocreales', 'Sordariomycetes',
     'Pezizomycotina', 'Ascomycota'),
    ('histoplasma_capsulatum', 'Histoplasma', 'Onygenales', 'Eurotiomycetes',
     'Pezizomycotina', 'Ascomycota'),
    ('kluyveromyces_lactis', 'Kluyveromyces', 'Saccharomycetales',
     'Saccharomycetes', 'Saccharomycotina', 'Ascomycota'),
    ('laccaria_bicolor', 'Laccaria', 'Agaricales', 'Agaricomycetes',
     'Agaricomycotina', 'Basidiomycota'),
    ('lodderomyces_elongisporus', 'Lodderomyces', 'Saccharomycetales',
     'Saccharomycetes', 'Saccharomycotina', 'Ascomycota'),
    ('magnaporthe_grisea', 'Magnaporthe', 'Magnaporthales',
     'Sordariomycetes', 'Pezizomycotina', 'Ascomycota'),
    ('neurospora_crassa', 'Neurospora', 'Sordariales', 'Sordariomycetes',
     'Pezizomycotina', 'Ascomycota'),
    ('phanerochaete_chrysosporium', 'Phanerochaete', 'Polyporales',
     'Agaricomycetes', 'Agaricomycotina', 'Basidiomycota'),
    ('pichia_stipitis', 'Pichia', 'Saccharomycetales', 'Saccharomycetes',
     'Saccharomycotina', 'Ascomycota'),
    ('rhizopus_oryzae', 'Rhizopus', 'Mucorales', '-', 'Mucoromycotina',
     'Mucoromycota'),
    ('saccharomyces_cerevisiae_S288C', 'Saccharomyces', 'Saccharomycetales',
     'Saccharomycetes', 'Saccharomycotina', 'Ascomycota'),
    ('saccharomyces_cerevisiae_rm11-1a_1', 'Saccharomyces', 'Saccharomycetales',
     'Saccharomycetes', 'Saccharomycotina', 'Ascomycota'),
    ('schizosaccharomyces_pombe', 'Schizosaccharomyces',
     'Schizosaccharomycetales', 'Schizosaccharomycetes', 'Taphrinomycotina',
     'Ascomycota'),
    ('ustilago_maydis', 'Ustilago', 'Ustilaginales', 'Ustilaginomycetes',
     'Ustilaginomycotina', 'Basidiomycota'),
    ('yarrowia_lipolytica', 'Yarrowia', 'Saccharomycetales', 'Saccharomycetes',
     'Saccharomycotina', 'Ascomycota')
]


def main():
    '''Main function'''
    argparse_usage = 'get_augustus_species.py -g <genus_name>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-g', '--genus_name', nargs=1, required=True,
        help='Genus name'
    )
    parser.add_argument(
        '-e', '--email_address', nargs=1, required=True,
        help='E-mail address for Entrez usage'
    )

    args = parser.parse_args()
    genus_name = args.genus_name[0]
    email_address = args.email_address[0]

    # Register E-mail address
    Entrez.email = email_address

    # Run functions :) Slow is as good as Fast
    get_augustus_species(genus_name)


def get_augustus_species(genus_name):
    '''Parse species data'''
    d_species = defaultdict(list)
    for tup in SPECIES_DATA:
        augustus_species, genus, order, tax_class, subphylum, phylum = tup
        for element in [genus, order, tax_class, subphylum, phylum]:
            d_species[element].append(augustus_species)

    # Get taxonomic info for input genus name
    print('Validating input genus name...')

    # Get taxonomy info from NCBI taxonomy
    taxon2 = '"' + genus_name + '"'
    handle = Entrez.esearch(
        db='taxonomy', term=taxon2, rettype='gb', retmode='text'
    )
    record = Entrez.read(handle, validate=False)
    handle.close()

    if record['IdList'] != []:
        handle2 = Entrez.efetch(
            db='taxonomy', id=record['IdList'][0], retmode='xml'
        )
        record2 = Entrez.read(handle2, validate=False)
        handle2.close()

    else:
        sys.exit(
            '[ERROR] The taxon "{}" you provided is invalid. Please check NCBI '
            'Taxonomy'.format(genus)
        )

    rank = record2[0]['Rank']
    lineage = record2[0]['Lineage']
    tax_list = record2[0]['LineageEx']

    # Initialization
    order = ''
    tax_class = ''
    subphylum = ''
    phylum = ''
    for tax_element in tax_list:
        if tax_element['Rank'] == 'phylum':
            phylum = tax_element['ScientificName']
        elif tax_element['Rank'] == 'subphylum':
            subphylum = tax_element['ScientificName']
        elif (
                tax_element['Rank'] == 'no rank' and
                re.search(r'.*cotina', tax_element['ScientificName'])
        ):
            subphylum = tax_element['ScientificName']
        elif tax_element['Rank'] == 'class':
            tax_class = tax_element['ScientificName']
        elif tax_element['Rank'] == 'order':
            order = tax_element['ScientificName']

    print('\n===')
    print('Taxon: %s' % (genus_name))
    print('Rank: %s' % (rank))
    print('Lineage: %s' % (lineage))
    print('===\n')

    print('Suggestions:')
    if genus_name in d_species:
        print('We found them in genus level: {}'.format(genus_name))
        print('\n'.join(d_species[genus_name]))
    elif order in d_species:
        print('We found them in order level: {}'.format(order))
        print('\n'.join(d_species[order]))
    elif subphylum in d_species:
        print('We found them in subphylum level: {}'.format(subphylum))
        print('\n'.join(d_species[subphylum]))
    elif phylum in d_species:
        print('We found them in phylum level: {}'.format(phylum))
        print('\n'.join(d_species[phylum]))
    else:
        print('We did not find any augustus species for the input "{}"'.format(
            genus_name
        ))


if __name__ == '__main__':
    main()
