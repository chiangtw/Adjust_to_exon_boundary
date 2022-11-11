#! /usr/bin/env python

import argparse
from circmimi.annotation import Annotation as AnnotationBase


class Annotation(AnnotationBase):
    def _get_nearest_junc_site(self, JuncSiteType, chr_name, pos, strand=None, dist=None):
        query_result = self.session\
            .query(JuncSiteType)\
            .filter_by(chr_id=self.chr_dict.get(chr_name))

        if strand is not None:
            query_result = query_result\
                .filter_by(strand_id=self.strand_dict.get(strand))

        up_result = query_result.filter(JuncSiteType.junc_site >= pos).first()
        down_result = query_result\
            .filter(JuncSiteType.junc_site <= pos)\
            .order_by(JuncSiteType.id.desc()).first()

        result = self._get_nearest_site(up_result, down_result, pos)

        if result and dist:
            if abs(result.junc_site - pos) > dist:
                result = None

        return result

    @staticmethod
    def _get_nearest_site(up_result, down_result, pos):
        if up_result and down_result:
            result = min(
                [
                    up_result,
                    down_result
                ],
                key=lambda JuncSite: abs(JuncSite.junc_site - pos)
            )

        elif not up_result and down_result:
            result = down_result

        elif up_result and not down_result:
            result = up_result

        else:
            result = None

        return result


def read_region_file(region_file_object):
    for line in region_file_object:
        data = line.rstrip('\n').split('\t')
        chrm, pos1, pos2, strand = data[:4]
        pos1 = int(pos1)
        pos2 = int(pos2)

        yield chrm, pos1, pos2, strand, data


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('anno_db')
    parser.add_argument(
        'region_file',
        type=argparse.FileType('r'),
        help="4-columns TSV file: (chrm, pos1, pos2, strand)"
    )
    parser.add_argument('--dist', type=int)
    parser.add_argument('--na_value', default='')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    anno_db = Annotation(args.anno_db)

    for chrm, pos1, pos2, strand, data in read_region_file(args.region_file):

        if strand == "+":
            donor_pos, acceptor_pos = pos2, pos1
        elif strand == "-":
            donor_pos, acceptor_pos = pos1, pos2

        nearest_donor = anno_db.get_nearest_donor_site(chrm, donor_pos, strand, args.dist + 1)
        nearest_acceptor = anno_db.get_nearest_acceptor_site(chrm, acceptor_pos, strand, args.dist + 1)

        if nearest_donor:
            nearest_donor_pos = nearest_donor.junc_site
            donor_adjust_bases = nearest_donor_pos - donor_pos
        else:
            nearest_donor_pos = args.na_value
            donor_adjust_bases = args.na_value

        if nearest_acceptor:
            nearest_acceptor_pos = nearest_acceptor.junc_site
            acceptor_adjust_bases = nearest_acceptor_pos - acceptor_pos
        else:
            nearest_acceptor_pos = args.na_value
            acceptor_adjust_bases = args.na_value

        if strand == "+":
            adjust_pos1 = nearest_acceptor_pos
            adjust_pos2 = nearest_donor_pos
            pos1_adjust_bases = acceptor_adjust_bases
            pos2_adjust_bases = donor_adjust_bases
        elif strand == "-":
            adjust_pos1 = nearest_donor_pos
            adjust_pos2 = nearest_acceptor_pos
            pos1_adjust_bases = donor_adjust_bases
            pos2_adjust_bases = acceptor_adjust_bases

        print(*data, adjust_pos1, adjust_pos2, pos1_adjust_bases, pos2_adjust_bases, sep='\t', flush=True)
