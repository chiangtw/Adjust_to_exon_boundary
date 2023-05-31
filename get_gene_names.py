#! /usr/bin/env python

import argparse
from circmimi.annotation import Annotation


def read_region_file(region_file_object):
    for line in region_file_object:
        data = line.rstrip('\n').split('\t')
        chrm, pos1, pos2, strand = data[:4]
        pos1 = int(pos1)
        pos2 = int(pos2)

        yield chrm, pos1, pos2, strand, data


def get_gene_names(chrm, pos1, pos2, strand, anno_db, na_value):
    if strand == "+":
        donor_pos, acceptor_pos = pos2, pos1
    elif strand == "-":
        donor_pos, acceptor_pos = pos1, pos2

    donor = anno_db.get_donor_site(chrm, donor_pos, strand)
    acceptor = anno_db.get_acceptor_site(chrm, acceptor_pos, strand)

    donor_gene_set = set()
    if donor:
        for exon in donor.exons:
            for transcript in exon.transcript:
                donor_gene_set.add(transcript.gene.gene_symbol)
        donor_genes = ','.join(sorted(donor_gene_set))
    else:
        donor_genes = na_value

    acceptor_gene_set = set()
    if acceptor:
        for exon in acceptor.exons:
            for transcript in exon.transcript:
                acceptor_gene_set.add(transcript.gene.gene_symbol)
        acceptor_genes = ','.join(sorted(acceptor_gene_set))
    else:
        acceptor_genes = na_value

    if (len(donor_gene_set) > 0) and (len(acceptor_gene_set) > 0):
        common_genes = donor_gene_set & acceptor_gene_set
        if len(common_genes) > 0:
            intragenic = 1
        else:
            intragenic = 0
    else:
        intragenic = na_value

    if strand == "+":
        pos1_genes = acceptor_genes
        pos2_genes = donor_genes
    elif strand == "-":
        pos1_genes = donor_genes
        pos2_genes = acceptor_genes

    return pos1_genes, pos2_genes, intragenic


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('anno_db')
    parser.add_argument(
        'region_file',
        type=argparse.FileType('r'),
        help="4-columns TSV file: (chrm, pos1, pos2, strand)"
    )
    parser.add_argument(
        '--na_value',
        default='NA',
        help='Placeholder for na values.'
    )

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    anno_db = Annotation(args.anno_db)

    for chrm, pos1, pos2, strand, data in read_region_file(args.region_file):
        pos1_genes, pos2_genes, intragenic = get_gene_names(chrm, pos1, pos2, strand, anno_db, args.na_value)
        print(*data, pos1_genes, pos2_genes, intragenic, sep='\t', flush=True)
