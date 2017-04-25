#!/usr/bin/python

'''
Download protein sequences of sister organisms from NCBI

Input: taxon name (e.g., Neurospora)
Output: a directory containing protein FASTA files gzipped
'''

# Import modules
import sys
import re
import os
from glob import glob
from Bio import Entrez
from random import sample
from argparse import ArgumentParser

# Parameters
ftp_base = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all'


# Main function
def main(argv):
    optparse_usage = (
        "downlaod_sister_orgs.py -d <download_dir> "
        "-t <taxon> -e <email_address>"
    )
    parser = ArgumentParser(usage=optparse_usage)
    parser.add_argument(
        "-d", "--download_dir", dest="download_dir", nargs=1,
        help="Download directory"
    )
    parser.add_argument(
        "-t", "--taxon", dest="taxon", nargs=1,
        help=(
            "Taxon that you want to download. You can choose any clade "
            "registered in NCBI (default: Fungi), but genus name is optimal "
            "(e.g., Neurospora)"
        )
    )
    parser.add_argument(
        "-n", "--num_sisters", dest="num_sisters", nargs=1,
        help="Number of sister orgasnisms"
    )
    parser.add_argument(
        "-e", "--email_address", dest="email_address", nargs=1,
        help="E-mail address for Entrez usage"
    )

    args = parser.parse_args()
    if args.download_dir:
        download_dir = os.path.abspath(args.download_dir[0])
    else:
        print '[ERROR] Please provide DOWNLOAD DIRECTORY'
        sys.exit(2)

    if args.taxon:
        taxon = args.taxon[0]
    else:
        print '[ERROR] Please provide TAXON'
        sys.exit(2)

    if args.num_sisters:
        num_sisters = int(args.num_sisters[0])
    else:
        print '[ERROR] Please provide NUBER OF SISTERS'
        sys.exit(2)

    if args.email_address:
        email_address = args.email_address[0]
    else:
        print '[ERROR] Please provide E-MAIL ADDRESS'
        sys.exit(2)

    # Register E-mail address
    Entrez.email = email_address

    # Run functions :) Slow is as good as Fast
    create_dir(download_dir)
    asm_ids = validate_taxon(taxon, num_sisters)
    download_genome(download_dir, taxon, asm_ids, num_sisters)


def create_dir(download_dir):
    if not os.path.exists(download_dir):
        os.mkdir(download_dir)


def validate_taxon(taxon, num_sisters):
    print 'Validate input taxon...'

    # Get taxonomy info from NCBI taxonomy
    taxon2 = '"' + taxon + '"'
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
            "Please check NCBI Taxonomy" % (taxon)
        )
        sys.exit(2)

    rank = record2[0]["Rank"]
    lineage = record2[0]["Lineage"]
    tax_list = record2[0]['LineageEx']

    # Initialization
    genus = ''
    family = ''
    order = ''
    Class = ''
    subphylum = ''
    phylum = ''
    kingdom = ''
    for tax_element in tax_list:
        if tax_element['Rank'] == 'kingdom':
            kingdom = tax_element['ScientificName']
        elif tax_element['Rank'] == 'phylum':
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
        elif tax_element['Rank'] == 'family':
            family = tax_element['ScientificName']
        elif tax_element['Rank'] == 'genus':
            genus = tax_element['ScientificName']

    print '\n==='
    print 'Taxon: %s' % (taxon)
    print 'Rank: %s' % (rank)
    print 'Lineage: %s' % (lineage)
    print '===\n'

    asm_ids = []
    i = 0
    taxa = [genus, family, order, Class, subphylum, phylum, kingdom]
    while len(asm_ids) < (num_sisters * 10):
        if not taxa[i]:
            i += 1
            continue

        # Get assembly IDs from NCBI
        handle3 = Entrez.esearch(
            db="assembly", term='%s[taxonomy]' % (taxa[i]),
            retmode="xml", retmax=1000000
        )
        record3 = Entrez.read(handle3, validate=False)
        handle3.close()

        tmp_asm_ids = record3["IdList"]
        if len(asm_ids) + len(tmp_asm_ids) > (num_sisters * 10):
            asm_ids_needed = (num_sisters * 10) - len(asm_ids)
            sampled_asm_ids = sample(tmp_asm_ids, asm_ids_needed)
            asm_ids.extend(sampled_asm_ids)

        else:
            asm_ids.extend(tmp_asm_ids)

        i += 1

    return asm_ids


def download_genome(download_dir, taxon, asm_ids, num_sisters):
    # Get Genbank accession ID and try to download .faa file
    cwd = os.getcwd()
    os.chdir(download_dir)
    print 'Downloading protein sequence files...'
    outtable = 'sister_orgs.list'
    outhandle = open(outtable, 'w')
    header_txt = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        'asm_id', 'organism', 'genbank_acc', 'kingdom', 'phylum',
        'subphylum', 'class', 'order', 'family', 'file_name'
    )
    outhandle.write(header_txt)
    out_faa = '-'
    num_downloaded_files = 0
    for asm_id in asm_ids:
        handle = Entrez.esummary(db='assembly', id=asm_id, retmode="xml")
        record = Entrez.read(handle, validate=False)
        handle.close()

        genbank_acc = (
            record['DocumentSummarySet']["DocumentSummary"][0]
            ["AssemblyAccession"]
        )

        asm_name = (
            record['DocumentSummarySet']["DocumentSummary"][0]["AssemblyName"]
        )
        asm_name2 = asm_name.replace(' ', '_')

        org_name = (
            record['DocumentSummarySet']['DocumentSummary'][0]['Organism']
        )

        tax_id = (
            record['DocumentSummarySet']['DocumentSummary'][0]['Taxid']
        )
        tax_tup = get_taxonomy(tax_id)
        kingdom, phylum, subphylum, Class, order, family, genus = tax_tup

        out_faa = os.path.join(
            download_dir, '%s_protein.faa.gz' % (genbank_acc)
        )
        out_faa_u = os.path.join(
            download_dir, '%s_protein.faa' % (genbank_acc)  # Unzipped file
        )
        if not os.path.exists(out_faa) and not os.path.exists(out_faa_u):
            acc_part1 = genbank_acc[0:3]
            acc_part2 = genbank_acc[4:7]
            acc_part3 = genbank_acc[7:10]
            acc_part4 = genbank_acc[10:13]
            url = '%s/%s/%s/%s/%s/%s_%s/*protein.faa.gz' % (
                ftp_base, acc_part1, acc_part2, acc_part3, acc_part4,
                genbank_acc, asm_name2
            )
            command = 'wget --quiet -nc %s' % (url)
            print '[Run] %s' % (command)
            os.system(command)
            downloaded_file = glob("%s_*protein.faa.gz" % (genbank_acc))
            if downloaded_file:
                command2 = 'mv %s %s' % (
                    downloaded_file[0], os.path.join(download_dir, out_faa)
                )
                print '[Run] %s' % (command2)
                os.system(command2)
                num_downloaded_files += 1
                if num_downloaded_files == num_sisters:
                    break
            else:
                print 'No downloadable file for this entry: %s' % (genbank_acc)
                out_faa = 'NA'

        # Write to table
        row_txt = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n' % (
            asm_id, org_name, genbank_acc, kingdom, phylum, subphylum, Class,
            order, family, out_faa
        )
        outhandle.write(row_txt)

    outhandle.close()
    os.chdir(cwd)
    print "\nDone. Check 'sister_orgs.list' for downloaded files details"


def get_taxonomy(tax_id):
    handle = Entrez.efetch(
        db="taxonomy", id=tax_id, retmode="xml"
    )
    record = Entrez.read(handle, validate=False)
    handle.close()

    # Initialization
    kingdom = ''
    phylum = ''
    subphylum = ''
    Class = ''
    order = ''
    family = ''
    genus = ''

    for rec in record[0]['LineageEx']:
        if rec['Rank'] == 'kingdom':
            kingdom = rec['ScientificName']
        elif rec['Rank'] == 'phylum':
            phylum = rec['ScientificName']
        elif rec['Rank'] == 'subphylum':
            subphylum = rec['ScientificName']
        elif (
            rec['Rank'] == 'no rank' and
            re.search(r'.*cotina', rec['ScientificName'])
        ):
            subphylum = rec['ScientificName']
        elif rec['Rank'] == 'class':
            Class = rec['ScientificName']
        elif rec['Rank'] == 'order':
            order = rec['ScientificName']
        elif rec['Rank'] == 'family':
            family = rec['ScientificName']
        elif rec['Rank'] == 'genus':
            genus = rec['ScientificName']

    tup = (kingdom, phylum, subphylum, Class, order, family, genus)

    return tup


if __name__ == "__main__":
    main(sys.argv[1:])
