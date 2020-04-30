#!/usr/bin/env python3

import sys

import argparse


def split_record(r, with_counts):
    r = r.split()
    if not with_counts:
        return [r[0].lower(), int(r[1]), int(r[2]), r[3].lower(), int(r[4]), int(r[5])]
    else:
        return [r[0].lower(), int(r[1]), int(r[2]), r[3].lower(), int(r[4]), int(r[5]), int(r[6])]


def check_records(list_of_records):
    chromosomes = set()
    for i in list_of_records:
        c1, a1, b1, c2, a2, b2, *_ = i
        chromosomes.add(c1)
        chromosomes.add(c2)
        if not a1 < b1 < a2 < b2:
            s = '\t'.join([str(j) for j in i])
            raise ValueError(f'Invalid coords in record: {s}')
    if len(list_of_records) == 0:
        raise ValueError("Empty list of records! (Maybe empty file?)")
    if len(chromosomes) > 1:
        raise ValueError('Interactions from different chromosomes!')


def convert(file_desc, res, with_counts, begin):
    records = [split_record(i, with_counts) for i in file_desc]
    records.sort(key=lambda x: x[1])
    check_records(records)
    bead_coords = set()
    for i in records:
        if not with_counts:
            _, a1, b1, _, a2, b2 = i
        else:
            _, a1, b1, _, a2, b2, c = i
        a = int(round(((a1 + b1) // 2 - begin) / res, 0))
        b = int(round(((a2 + b2) // 2 - begin) / res, 0))
        if b - a >= 2:
            if not with_counts:
                bead_coords.add((a, b))
            else:
                bead_coords.add((a, b, c))
    bead_coords = list(bead_coords)
    bead_coords.sort(key=lambda x: (x[0], x[1]))
    if with_counts:
        i = 0
        j = 1
        clean_bead_coords = []
        while j < len(bead_coords):
            a1, b1, c1 = bead_coords[i]
            a2, b2, c2 = bead_coords[j]
            if a1 == a2 and b1 == b2:
                while a1 == a2 and b1 == b2:
                    c1 += c2
                    j += 1
                    a2, b2, c2 = bead_coords[j]
                else:
                    clean_bead_coords.append((a1, b1, c1))
                    i = j
                    j = i + 1
            else:
                clean_bead_coords.append((a1, b1, c1))
                i += 1
                j = i + 1
        a1, b1, c1 = bead_coords[i]
        clean_bead_coords.append((a1, b1, c1))
        bead_coords = clean_bead_coords
    return bead_coords


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="translate genomic coords into model coords")
    parser.add_argument("resolution", type=int, help="Model resolution")
    parser.add_argument('begin', type=int, help='beginning of region of interest')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--with_counts', action='store_true', default=False, help="Include counts in 3rd column?")
    args = parser.parse_args()
    w = ''
    for i in convert(args.infile, args.resolution, args.with_counts, args.begin):
        if not args.with_counts:
            w += f':{i[0]}\t:{i[1]}\n'
        else:
            w += f':{i[0]}\t:{i[1]}\t{i[2]}\n'
    w = w[:-1]
    args.outfile.write(w)
