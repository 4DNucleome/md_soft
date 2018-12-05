import sys

import argparse


def split_record(r):
    r = r.split()
    return [r[0].lower(), int(r[1]), int(r[2]), r[3].lower(), int(r[4]), int(r[5])]


def check_records(list_of_records):
    chromosomes = set()
    for i in list_of_records:
        c1, a1, b1, c2, a2, b2 = i
        chromosomes.add(c1)
        chromosomes.add(c2)
        if not a1 < b1 < a2 < b2:
            s = '\t'.join([str(j) for j in i])
            raise ValueError(f'Invalid coords in record: {s}')
    if len(chromosomes) != 1:
        raise ValueError('Interactions from different chromosomes!')


def convert(file_desc, res):
    records = [split_record(i) for i in file_desc]
    records.sort(key=lambda x: x[1])
    check_records(records)
    begin = records[0][1]
    bead_coords = set()
    for i in records:
        _, a1, b1, _, a2, b2 = i
        a = int((a1 + (b1 - a1) / 2 - begin) / res) + 1
        b = int((a2 + (b2 - a2) / 2 - begin) / res) + 1
        if b - a >= 2:
            bead_coords.add((a, b))
    bead_coords = list(bead_coords)
    bead_coords.sort(key=lambda x: (x[0], x[1]))
    return bead_coords


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="translate genomic coords into model coords")
    parser.add_argument("resolution", type=int, help="Model resolution")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    w = ''
    for i in convert(args.infile, args.resolution):
        w += f'{i[0]}\t{i[1]}\n'
    w = w[:-1]
    args.outfile.write(w)
