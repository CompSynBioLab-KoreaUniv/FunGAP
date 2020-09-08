#!/usr/bin/env python3

'''
Download protein sequences of sister organisms from NCBI

Input: taxon name (e.g., Neurospora)
Output: a directory containing protein FASTA files gzipped
Last updated: Jul 13, 2020
'''

import os
import re
import sys
from argparse import ArgumentParser
from glob import glob
from random import sample

from Bio import Entrez

# Parameters
FTP_BASE = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all'


# Main function
def main():
    '''Main functin'''
    optparse_usage = (
        'downlaod_sister_orgs.py -d <download_dir> -t <taxon> '
        '-e <email_address>'
    )
    parser = ArgumentParser(usage=optparse_usage)
    parser.add_argument(
        '-d', '--download_dir', nargs='?', default='sister_orgs',
        help='Download directory'
    )
    parser.add_argument(
        '-t', '--taxon', nargs=1, required=True,
        help=(
            'Taxon that you want to download. You can choose any clade '
            'registered in NCBI, but genus name is optimal (e.g., Neurospora)'
        )
    )
    parser.add_argument(
        '-n', '--num_sisters', nargs='?', default=3, type=int,
        help='Number of sister organisms (default: 3)'
    )
    parser.add_argument(
        '-e', '--email_address', nargs=1, required=True,
        help='E-mail address for Entrez usage'
    )

    args = parser.parse_args()
    download_dir = os.path.abspath(args.download_dir)
    taxon = args.taxon[0]
    num_sisters = args.num_sisters
    email_address = args.email_address[0]

    # Register E-mail address
    Entrez.email = email_address

    # Run functions :) Slow is as good as Fast
    create_dir(download_dir)
    asm_ids = validate_taxon(taxon, num_sisters)
    download_genome(download_dir, asm_ids, num_sisters)


def create_dir(download_dir):
    '''Create directory'''
    if not os.path.exists(download_dir):
        os.mkdir(download_dir)


def validate_taxon(taxon, num_sisters):
    '''Validate taxon'''
    print('Validate input taxon...')

    # Get taxonomy info from NCBI taxonomy
    taxon2 = '"' + taxon + '"'
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
            '[ERROR] The taxon "{}" you provided is invalid. '
            'Please check NCBI Taxonomy'.format(taxon)
        )

    rank = record2[0]['Rank']
    lineage = record2[0]['Lineage']
    tax_list = record2[0]['LineageEx']

    # Initialization
    genus = ''
    family = ''
    order = ''
    tax_class = ''
    subphylum = ''
    phylum = ''
    kingdom = ''
    for tax_element in tax_list:
        cond1 = tax_element['Rank'] == 'no rank'
        cond2 = re.search(r'.*cotina', tax_element['ScientificName'])
        if tax_element['Rank'] == 'kingdom':
            kingdom = tax_element['ScientificName']
        elif tax_element['Rank'] == 'phylum':
            phylum = tax_element['ScientificName']
        elif tax_element['Rank'] == 'subphylum':
            subphylum = tax_element['ScientificName']
        elif cond1 and cond2:
            subphylum = tax_element['ScientificName']
        elif tax_element['Rank'] == 'class':
            tax_class = tax_element['ScientificName']
        elif tax_element['Rank'] == 'order':
            order = tax_element['ScientificName']
        elif tax_element['Rank'] == 'family':
            family = tax_element['ScientificName']
        elif tax_element['Rank'] == 'genus':
            genus = tax_element['ScientificName']

    print('\n===')
    print('Taxon: {}'.format(taxon))
    print('Rank: {}'.format(rank))
    print('Lineage: {}'.format(lineage))
    print('===\n')

    asm_ids = []
    i = 0
    taxa = [genus, family, order, tax_class, subphylum, phylum, kingdom]
    while len(asm_ids) < (num_sisters * 10):
        if not taxa[i]:
            i += 1
            continue

        # Get assembly IDs from NCBI
        handle3 = Entrez.esearch(
            db='assembly', term='{}[taxonomy]'.format(taxa[i]),
            retmode='xml', retmax=1000000
        )
        record3 = Entrez.read(handle3, validate=False)
        handle3.close()

        tmp_asm_ids = record3['IdList']
        if len(asm_ids) + len(tmp_asm_ids) > (num_sisters * 10):
            asm_ids_needed = (num_sisters * 10) - len(asm_ids)
            sampled_asm_ids = sample(tmp_asm_ids, asm_ids_needed)
            asm_ids.extend(sampled_asm_ids)

        else:
            asm_ids.extend(tmp_asm_ids)

        i += 1
    return asm_ids


def download_genome(download_dir, asm_ids, num_sisters):
    '''Get Genbank accession ID and try to download .faa file'''
    cwd = os.getcwd()
    os.chdir(download_dir)
    print('Downloading protein sequence files...')
    outtable = 'sister_orgs.list'
    outhandle = open(outtable, 'w')
    header_txt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        'asm_id', 'organism', 'genbank_acc', 'kingdom', 'phylum',
        'subphylum', 'class', 'order', 'family', 'file_name'
    )
    outhandle.write(header_txt)
    out_faa = '-'
    num_downloaded_files = 0
    for asm_id in asm_ids:
        handle = Entrez.esummary(db='assembly', id=asm_id, retmode='xml')
        record = Entrez.read(handle, validate=False)
        handle.close()

        genbank_acc = (
            record['DocumentSummarySet']['DocumentSummary'][0]
            ['AssemblyAccession']
        )

        asm_name = (
            record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']
        )
        asm_name2 = asm_name.replace(' ', '_')

        org_name = (
            record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
        )

        tax_id = (
            record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']
        )
        tax_tup = get_taxonomy(tax_id)
        kingdom, phylum, subphylum, tax_class, order, family = tax_tup[:6]

        out_faa = os.path.join(
            download_dir, '{}_protein.faa.gz'.format(genbank_acc)
        )
        out_faa_u = os.path.join(
            download_dir, '{}_protein.faa'.format(genbank_acc)  # Unzipped file
        )
        if not os.path.exists(out_faa) and not os.path.exists(out_faa_u):
            acc_part1 = genbank_acc[0:3]
            acc_part2 = genbank_acc[4:7]
            acc_part3 = genbank_acc[7:10]
            acc_part4 = genbank_acc[10:13]
            url = '{}/{}/{}/{}/{}/{}_{}/*protein.faa.gz'.format(
                FTP_BASE, acc_part1, acc_part2, acc_part3, acc_part4,
                genbank_acc, asm_name2
            )
            command = 'wget --quiet -nc {}'.format(url)
            print('[Run] {}'.format(command))
            os.system(command)
            downloaded_file = glob('{}_*protein.faa.gz'.format(genbank_acc))
            if downloaded_file:
                command2 = 'mv {} {}'.format(
                    downloaded_file[0], os.path.join(download_dir, out_faa)
                )
                print('[Run] {}'.format(command2))
                os.system(command2)
                num_downloaded_files += 1
                if num_downloaded_files == num_sisters:
                    break
            else:
                print('No downloadable file for this entry: {}'.format(
                    genbank_acc
                ))
                out_faa = 'NA'

        # Write to table
        row_txt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            asm_id, org_name, genbank_acc, kingdom, phylum, subphylum,
            tax_class, order, family, out_faa
        )
        outhandle.write(row_txt)

    outhandle.close()
    os.chdir(cwd)
    print('\nDone. Check "sister_orgs.list" for downloaded files details')


def get_taxonomy(tax_id):
    '''Get taxonomy'''
    handle = Entrez.efetch(db='taxonomy', id=tax_id, retmode='xml')
    record = Entrez.read(handle, validate=False)
    handle.close()

    # Initialization
    kingdom, phylum, subphylum, tax_class, order, family, genus = [''] * 7
    for rec in record[0]['LineageEx']:
        cond1 = rec['Rank'] == 'no rank'
        cond2 = re.search(r'.*cotina', rec['ScientificName'])
        if rec['Rank'] == 'kingdom':
            kingdom = rec['ScientificName']
        elif rec['Rank'] == 'phylum':
            phylum = rec['ScientificName']
        elif rec['Rank'] == 'subphylum':
            subphylum = rec['ScientificName']
        elif cond1 and cond2:
            subphylum = rec['ScientificName']
        elif rec['Rank'] == 'class':
            tax_class = rec['ScientificName']
        elif rec['Rank'] == 'order':
            order = rec['ScientificName']
        elif rec['Rank'] == 'family':
            family = rec['ScientificName']
        elif rec['Rank'] == 'genus':
            genus = rec['ScientificName']

    tup = (kingdom, phylum, subphylum, tax_class, order, family, genus)
    return tup


if __name__ == '__main__':
    main()
