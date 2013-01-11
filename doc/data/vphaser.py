#!/usr/bin/env python3

import sys

from os import close, remove
from os.path import exists, join
from argparse import ArgumentParser, FileType
from subprocess import Popen
from tempfile import mkstemp
from textwrap import dedent

from Bio import SeqIO


PERL = '/usr/bin/perl'
RC454_DIR = '/home/sci/uds/RC454_SoftwarePackage'
VPHASER_DIR = '/home/sci/uds/VpSoftwarePackage'


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = ArgumentParser(description='run RC454')
    parser.add_argument('--perl', default=PERL)
    parser.add_argument('--rc454', default=RC454_DIR)
    parser.add_argument('--vphaser', default=VPHASER_DIR)
    parser.add_argument('reads', type=FileType('r'))
    parser.add_argument('quals', type=FileType('r'))
    parser.add_argument('assembly', type=FileType('r'))
    parser.add_argument('--primers')
    parser.add_argument('--genelist')
    parser.add_argument('--noendvariant', type=int, default=10)
    subpar = parser.add_mutually_exclusive_group()
    subpar.add_argument('--regions')
    subpar.add_argument('--autoregions', action='store_true', default=False)
    parser.add_argument('output')

    ns = parser.parse_args(args)

    runmosaik = join(ns.rc454, 'runMosaik.pl')
    rc454 = join(ns.rc454, 'rc454.pl')
    vphaser = join(ns.vphaser, 'vphaser.pl')
    vprofiler = join(ns.vphaser, 'vprofiler.pl')

    # create temporary files so that the variables are guaranteed to be around for the finally: block
    fd, genelist = mkstemp(); close(fd)
    fd, prefix = mkstemp(); close(fd)
    qlxfile = prefix + '.qlx'
    fd, vprofile = mkstemp(); close(fd)

    try:
        record = SeqIO.read(ns.assembly, 'fasta')

        if ns.genelist is None:
            with open(genelist, 'w') as fh:
                fh.write('%s\t%s\t%s\n' % (record.id, 1, len(record)))
        else:
            genelist = ns.genelist

        if ns.autoregions:
            regions = genelist
        else:
            regions = ns.regions

        # perl ../runMosaik.pl VTest_reads.fa VTest_reads.qual VTest_assembly.fa VTest_raw
        p = Popen(
            [
                ns.perl,
                runmosaik,
                ns.reads.name,
                ns.quals.name,
                ns.assembly.name,
                prefix,
                '-qlxonly'
            ],
            stderr=sys.stderr,
            stdout=sys.stdout
        )
        p.wait()

        # perl ../rc454.pl VTest_raw.qlx VTest_reads.fa VTest_reads.qual VTest_assembly.fa rcVTest -bam -primers=VTest_primers.txt -genelist=VTest_genelist.txt
        args = [
            ns.perl,
            rc454,
            qlxfile,
            ns.reads.name,
            ns.quals.name,
            ns.assembly.name,
            ns.output,
            # '-bam',
            '-genelist=%s' % genelist
        ]

        # handle the primers case
        if ns.primers is not None:
            args.insert(-1, '-primers=%s' % ns.primers)

        p = Popen(
            args,
            stderr=sys.stderr,
            stdout=sys.stdout
        )
        p.wait()

        finalqlxfile = ns.output + '_final.qlx'

        # perl ../vphaser.pl -i rcVTest_final.qlx -o vp_VTest
        p = Popen(
            [
                ns.perl,
                vphaser,
                '-i', finalqlxfile,
                '-o', ns.output
            ],
            stderr=sys.stderr,
            stdout=sys.stdout
        )
        p.wait()

        with open(vprofile, 'w') as fh:
            fh.write(dedent('''\
                >Qlx
                %s\t%s

                >VPhaser
                %s\t%s

                >Consensus
                %s\t%s

                >Genelist
                %s\t%s
                ''' % (
                    finalqlxfile, record.id,
                    ns.output + '_calls.txt', record.id,
                    ns.assembly.name, record.id,
                    genelist, record.id
                )
            ))

        args = [
            ns.perl,
            vprofiler,
            '-i', vprofile,
            '-o', ns.output,
            '-noendvariant=%d' % ns.noendvariant,
            '-nt',
            '-codon'
        ]

        if regions is not None:
            with open(vprofile, 'a') as fh:
                fh.write('\n>Region\n%s\t%s' % (regions, record.id))
            args.extend(['-haplo', '-haploseq'])

        # perl ../vprofiler.pl -i vprofiler_input_VTest.txt -o vpro -noendvariant=10 -nt -codon -haplo -haploseq
        p = Popen(
            args,
            stderr=sys.stderr,
            stdout=sys.stdout
        )
        p.wait()

    finally:
        if ns.genelist is None and exists(genelist):
            remove(genelist)
        if exists(prefix):
            remove(prefix)
        if exists(qlxfile):
            remove(qlxfile)
        if exists(vprofile):
            remove(vprofile)

    return 0


if __name__ == '__main__':
    sys.exit(main())
