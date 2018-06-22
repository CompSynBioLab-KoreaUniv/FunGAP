#!/usr/bin/python

'''
Get Augustus species argument by taxon name
Author Byoungnam Min on May 29, 2018
'''

# Import modules
import re
import sys
from Bio import Entrez
from collections import defaultdict
from argparse import ArgumentParser

# Parameter
species_data = [
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
    ('candida_guilliermondii', 'Candida', 'Saccharomycetales', 'Saccharomycetes',
     'Saccharomycotina', 'Ascomycota'),
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


# Main function
def main(argv):
    argparse_usage = 'get_augustus_species.py -g <genus_name>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-g", "--genus_name", dest="genus_name", nargs=1,
        help="Genus name"
    )
    parser.add_argument(
        "-e", "--email_address", dest="email_address", nargs=1,
        help="E-mail address for Entrez usage"
    )

    args = parser.parse_args()
    if args.genus_name:
        genus_name = args.genus_name[0]
    else:
        print '[ERROR] Please provide INPUT FASTA'
        parser.print_help()
        sys.exit(2)

    if args.email_address:
        email_address = args.email_address[0]
    else:
        print '[ERROR] Please provide E-MAIL ADDRESS'
        sys.exit(2)

    # Register E-mail address
    Entrez.email = email_address

    # Run functions :) Slow is as good as Fast
    get_augustus_species(genus_name)


def get_augustus_species(genus_name):
    # Parse species data
    D_species = defaultdict(list)
    for tup in species_data:
        augustus_species, genus, order, Class, subphylum, phylum = tup
        for element in [genus, order, Class, subphylum, phylum]:
            D_species[element].append(augustus_species)

    # Get taxonomic info for input genus name
    print 'Validating input genus name...'

    # Get taxonomy info from NCBI taxonomy
    taxon2 = '"' + genus_name + '"'
    handle = Entrez.esearch(
        db="taxonomy", term=taxon2, rettype="gb", retmode="text"
    )
    record = Entrez.read(handle, validate=False)
    handle.close()

    if record['IdList'] != []:
        handle2 = Entrez.efetch(
            db="taxonomy", id=record['IdList'][0], retmode="xml"
        )
        record2 = Entrez.read(handle2, validate=False)
        handle2.close()

    else:
        print (
                "The taxon '%s' you provided is invalid. "
                "Please check NCBI Taxonomy" % (genus)
        )
        sys.exit(2)

    rank = record2[0]["Rank"]
    lineage = record2[0]["Lineage"]
    tax_list = record2[0]['LineageEx']

    # Initialization
    order = ''
    Class = ''
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
            Class = tax_element['ScientificName']
        elif tax_element['Rank'] == 'order':
            order = tax_element['ScientificName']

    print '\n==='
    print 'Taxon: %s' % (genus_name)
    print 'Rank: %s' % (rank)
    print 'Lineage: %s' % (lineage)
    print '===\n'

    print 'Suggestions:'
    if genus_name in D_species:
        print 'We found them in genus level: {}'.format(genus_name)
        print '\n'.join(D_species[genus_name])
    elif order in D_species:
        print 'We found them in order level: {}'.format(order)
        print '\n'.join(D_species[order])
    elif subphylum in D_species:
        print 'We found them in subphylum level: {}'.format(subphylum)
        print '\n'.join(D_species[subphylum])
    elif phylum in D_species:
        print 'We found them in phylum level: {}'.format(phylum)
        print '\n'.join(D_species[phylum])
    else:
        print 'We did not find any augustus species for the input "{}"'.format(
            genus_name
        )


if __name__ == "__main__":
    main(sys.argv[1:])