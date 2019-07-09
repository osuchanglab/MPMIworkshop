#!/usr/bin/env python2.7
from reportlab.lib import colors
# from reportlab.lib.units import cm
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import inch
from reportlab.lib.colors import HexColor
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio.SeqFeature import SeqFeature, FeatureLocation
# from math import log1p as log
# from math import ceil
from collections import defaultdict
# import sys
import csv
import argparse
import logging
import os

def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def init_logger(args):
    logger = logging.getLogger(__name__)
    ch = logging.StreamHandler()
    if args.debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
        ch.setLevel(logging.WARNING)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return(logger)


parser = argparse.ArgumentParser(
        description='Generate map image from a genbank file and feature file.')
parser.add_argument(
        'gbk', help='Multi-Genbank file with genomes of interest.',
        type=extant_file)
parser.add_argument(
        'featurefile', help='Features and coordinates to include.',
        type=extant_file)
parser.add_argument(
        '--colorfile', help='Custom colors to use for each feature group.',
        type=extant_file)
parser.add_argument(
        '--color', help='Default color for features. (default = gray)',
        type=str, default='gray')
parser.add_argument(
        '--blast', help='File with BLAST search results for connecting hits.',
        type=extant_file)
parser.add_argument(
        '--outname', help='Name out output file to print (default = '
                          'input.gbk.pdf)', type=str)
parser.add_argument(
        '--outtype', help='Type of output file to print (default = PDF)',
        choices=['PDF', 'SVG', 'EPS'], default='PDF', type=str)
parser.add_argument(
        '--size', help='Size of output.',
        choices=['mbio', 'mbio_half', 'mbio_half_height',
                 'mbio_half_width'],
        default='big', type=str)
parser.add_argument(
        '--reverse', help='File with contigs to reverse in output.',
        type=str)
parser.add_argument(
        '--hfactor', help='Modify height by this factor. default = 1.0',
        type=float, default=1.0)
parser.add_argument(
        '--wfactor', help='Modify width by this factor. default = 1.0',
        type=float, default=1.0)
parser.add_argument(
        '--orientation', help='Orientation for final output.',
        choices=['portrait', 'landscape'])
parser.add_argument(
        '--verbose', help='Print progress messages.', action='store_true')
parser.add_argument(
        '--debug', help='Print debugging messages.', action='store_true')
args = parser.parse_args()
args.logger = init_logger(args)
offset = {}
custom = {}

reverse = set()

if args.reverse:
    with open(args.reverse, 'rb') as revfh:
        for line in revfh:
            line = line.strip()
            reverse.add(line)

blastresults = defaultdict(list)
if args.blast:
    msg = 'Getting BLAST alignment data from {}.'
    args.logger.debug(msg.format(args.blast))
    with open(args.blast, 'rb') as tsv:
        for line in csv.reader(tsv, delimiter='\t'):
            # query id, subject id, % identity, alignment length, mismatches,
            # gap opens, q. start, q. end, s. start, s. end, evalue, bit score
            if line[0].startswith('#'):
                continue
            query = line[0]
            target = line[5]
            strand = line[4]
            pident = float(line[9])/float(line[10])*100
            qstart = int(line[2])
            qend = int(line[3])
            if strand == "-":
                temp = qstart
                qstart = qend
                qend = temp
            sstart = int(line[7])
            send = int(line[8])
            data = (target, pident, qstart, qend, sstart, send)
            blastresults[query].append(data)

featcolors = {}
if args.colorfile:
    msg = 'Getting color data from {}.'
    args.logger.debug(msg.format(args.colorfile))
    with open(args.colorfile, 'rb') as cf:
        for line in csv.reader(cf, delimiter='\t'):
            msg = 'Identified feature group:{} -> color:{}.'
            args.logger.debug(msg.format(line[0], line[1]))
            featcolors[line[0]] = line[1]

features = {}
coords = defaultdict(list)
with open(args.featurefile, 'rb') as ff:
    for line in csv.reader(ff, delimiter='\t'):
        # rec_name feat_type feat_group locus_tag/coordinates
        if line[1] == 'feature':
            msg = 'Feature:{} -> Group:{}'
            args.logger.debug(msg.format(line[3], line[2]))
            features[line[3]] = line[2]
        elif line[1] == 'coordinates':
            msg = 'Coordinates:{} -> Group:{}'
            args.logger.debug(msg.format(line[3], line[2]))
            coords[line[0]].append((line[2], line[3]))


gbk_rec = SeqIO.parse(args.gbk, 'gb')

name = os.path.basename(args.gbk)

gd_diagram = GenomeDiagram.Diagram(name)
max_len = 0

gbk_rec_to_parse = []

crosslinks = {}
tn = 1
space = 1
tracks = {}
for record in gbk_rec:
    msg = 'Working on genbank record {}'
    if record.name in reverse:
        print("yes its in reverse")
        record = record.reverse_complement(id=True,name=True)
    args.logger.info(msg.format(record.name))
    start = 0
    end = int(len(record.seq))
    max_len = max(max_len, end)
    gt_fontsize = 12
    alen = 0.33
    aheight = 0.33
    gd_track_for_features = gd_diagram.new_track(
            tn,
            name=record.name,
            greytrack=True,
            greytrack_fontsize=gt_fontsize,
            greytrack_labels=2,
            greytrack_font_colour=colors.black,
            greytrack_font_color=colors.black,
            scale_smalltick_interval=100000,
            scale=1,
            start=start,
            end=end,
            height=0.5,
            # track_size=0.75,
            level=1)
    gd_feature_set = gd_track_for_features.new_set()

    tracks[record.name] = tn

    coordlist = coords.get(record.name)
    if coordlist is not None:
        for coord in coords.get(record.name):
            ftype = coord[0]
            delimiter = ','
            if ',' in coord[1]:
                pass
            elif ':' in coord[1]:
                delimiter = ':'
            elif ';' in coord[1]:
                delimiter = ';'
            elif '..' in coord[1]:
                delimiter = '..'
            try:
                start, end = map(int, coord[1].split(delimiter))
            except:
                msg = 'Unable to split coordinates {} for record {}.'
                args.logger.warning(msg.format(coord[1], record.name))
                msg = 'These positions will be disregarded.'
                args.logger.warning(msg)
                continue
            try:
                fcolor = featcolors[ftype]
            except KeyError:
                msg = 'No color provided for type {}.'
                args.logger.debug(msg.format(ftype))
                fcolor = args.color
            if record.name in reverse:
                start = int(len(record.seq)) - start
                end = int(len(record.seq)) - end
                postemp = start
                start = end
                end = postemp
            msg = 'Adding {} colored coordinates box at {},{} for {}'
            args.logger.debug(msg.format(fcolor, start, end, record.name))
            if ftype == "IS":
                feat = SeqFeature(
                        FeatureLocation(start, end), strand=-1)
            else:
                feat = SeqFeature(
                        FeatureLocation(start, end), strand=1)
            gd_feature_set.add_feature(
                    feat, sigil='BOX', color=fcolor, label=False)

    for feature in record.features:
        if feature.type != 'CDS':
            # Exclude this feature
            continue
        if 'locus_tag' in feature.qualifiers:
            flocus = feature.qualifiers['locus_tag'][0]
        else:
            continue
        try:
            ftype = features[flocus]
        except KeyError:
            # msg = 'Locus {} not included.'
            # args.logger.debug(msg.format(flocus))
            continue
        try:
            fcolor = featcolors[ftype]
        except KeyError:
            msg = 'No color provided for type {}.'
            args.logger.debug(msg.format(ftype))
            fcolor = args.color
        msg = ('Adding {} colored feature box on strand {} for {} of type {} '
               'from {}.')
        args.logger.debug(
                msg.format(fcolor, feature.strand, flocus, ftype, record.name))
        gd_feature_set.add_feature(
                feature, sigil='BOX', color=fcolor, label=False)

    tn += space


revspace = space - 1
pidcolor = colors.red
color_range_min = 0
color_range_max = 100

for query in blastresults:
    msg = 'Found BLAST searches by query {}.'
    args.logger.debug(msg.format(query))
    for hit in blastresults[query]:
        target = hit[0]
        pident = hit[1]
        qstart = hit[2]
        qend = hit[3]
        sstart = hit[4]
        send = hit[5]
        if query in reverse:
            qstart = gd_diagram.tracks[tracks[query]].end - qstart
            qend = gd_diagram.tracks[tracks[query]].end - qend
        if target in reverse:
            sstart = gd_diagram.tracks[tracks[target]].end - sstart
            send = gd_diagram.tracks[tracks[target]].end - send
        msg = 'link:{}-{}_id{}_q{}:{}_s{}:{}'
        args.logger.debug(msg.format(query, target, pident, qstart, qend,
                          sstart, send))
        t1 = gd_diagram.tracks[tracks[query]]
        t2 = gd_diagram.tracks[tracks[target]]
        if pident > 99.9:
            color = HexColor("#910202")
        elif pident > 98:
            color = HexColor("#ff0000")
        elif pident > 90:
            color = HexColor("#fb4444")
        elif pident >80:
            color = HexColor("#ff8686")
        #elif pident > 50:
        else:
            color = HexColor("#fdc7c7")
        #else:
            #color = HexColor("#ffe2e2")
            #color = HexColor("#ffffff")
        link = CrossLink((t1, qstart, qend), (t2, sstart, send),
                         color)
                         #color, HexColor("#ff8686"))
        gd_diagram.cross_track_links.append(link)

gd_diagram.renumber_tracks(low=(tn-space), step=-space)

# orientation = 'landscape'
# if fragments > 1 or len(gd_diagram.tracks) > 1:
#     orientation = 'portrait'
# pagesize = (10*inch, 0.5*inch*fragments)
hfactor = 1.25*args.hfactor
wfactor = 1.25*args.wfactor
pagesize = (8.5*inch*wfactor, 11*inch*hfactor)
orientation = 'portrait'
mbioh = 9.0625*inch*hfactor
mbiow = 6.875*inch*wfactor

testing = 1
if testing:
    pagesize = (mbiow/2, mbioh/4)
    orientation = 'landscape'

if args.size == 'mbio':
    pagesize = (mbiow, mbioh)
    orientation = 'portrait'
elif args.size == 'mbio_half':
    pagesize = (mbiow/2, mbioh/2)
    orientation = 'portrait'
elif args.size == 'mbio_half_height':
    pagesize = (mbiow, mbioh/2)
    orientation = 'landscape'
elif args.size == 'mbio_half_width':
    pagesize = (mbiow/2, mbioh)
    orientation = 'portrait'
if args.orientation:
    orientation = args.orientation

# print 'fragments : {}'.format(fragments)
outname = '{}.{}'.format(name, args.outtype)

if args.outname:
    outname = args.outname

print 'Printing {} file to {}'.format(args.outtype, outname)
gd_diagram.draw(
        format='linear', pagesize=pagesize, orientation=orientation,
        fragments=1, start=0, end=float(max_len))
gd_diagram.write(outname, args.outtype)
