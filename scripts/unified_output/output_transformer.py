import argparse
import os
import sys
from unified.Event import Event
from gtf_utils.GTFParser import GTFParser
from tool_parser.ASGALParser import ASGALParser
from tool_parser.IRFinderParser import IRFinderParser
from tool_parser.MAJIQParser import MAJIQParser
from tool_parser.SpladderParser import SplAdderParser
from tool_parser.WhippetParser import WhippetParser
from tool_parser.ASpliParser import ASpliParser
from tool_parser.EventPointerParser import EventPointerParser
from tool_parser.SGSeqReader import SGSeqReader
from tool_parser.SGSeqAnnoReader import SGSeqAnnoReader
from unified.EventAnnotationReader import EventAnnotationReader
from unified.compare_files import compare_with_annotation


def run_compare(args):
    gtf = GTFParser(args.gtf)
    event_annotation = EventAnnotationReader(args.event_annotation, gtf, args.combine_me)

    compare_with_annotation(event_annotation, args.compare_file, args.stats_outfile, args.strict, args.threshold)


def run_create(args):
    if args.spladder_dir is not None:
        spladder = SplAdderParser(args.spladder_dir, args.combine_me)
        run(spladder, args)

    if args.majiq_dir is not None:
        psi_filepath = None
        voila_filepath = None
        for file in os.listdir(args.majiq_dir):
            if file.endswith("psi.tsv"):
                psi_filepath = os.path.join(args.majiq_dir, file)
            if file.endswith("voila.tsv"):
                voila_filepath = os.path.join(args.majiq_dir, file)
        if psi_filepath is None or voila_filepath is None:
            print("Did not find correct MAJIQ files; check that they are named correctly!")
            exit(1)

        majiq = MAJIQParser(psi_filepath, voila_filepath, args.outdir, args.combine_me)
        run(majiq, args)

    if args.whippet_file is not None:
        whippet = WhippetParser(args.whippet_file, args.combine_me)
        run(whippet, args)

    if args.asgal_file is not None:
        gtf = GTFParser(args.gtf)
        asgal = ASGALParser(args.asgal_file, gtf, args.combine_me)
        run(asgal, args)

    if args.irfinder_file is not None:
        irfinder = IRFinderParser(args.irfinder_file)
        run(irfinder, args)

    if args.aspli_dir is not None:
        gtf = GTFParser(args.gtf)
        discovery_file = os.path.join(args.aspli_dir, "as_discovery.tab")
        exon_counts_file = os.path.join(args.aspli_dir, "exon.counts.tab")
        intron_counts_file = os.path.join(args.aspli_dir, "intron.counts.tab")
        aspli = ASpliParser(discovery_file, exon_counts_file, intron_counts_file, gtf)
        run(aspli, args)

    if args.eventpointer_file is not None:
        gtf = GTFParser(args.gtf)
        eventpointer = EventPointerParser(args.eventpointer_file, gtf, args.combine_me)
        run(eventpointer, args)

    if args.sgseq_denovo is not None:
        gtf = GTFParser(args.gtf)
        sgseq_denovo = SGSeqReader(args.sgseq_denovo, gtf, args.combine_me)
        run(sgseq_denovo, args)

    if args.sgseq_anno is not None:
        gtf = GTFParser(args.gtf)
        sgseq_anno = SGSeqAnnoReader(args.sgseq_anno, gtf, args.combine_me)
        run(sgseq_anno, args)


def run(tool, args):
    header = "chr\tgene\tid\tstrand\tevent_type\tcount\tstart_coordinates\tend_coordinates\n"
    outfile = args.outdir + "/" + os.path.basename(os.path.normpath(args.outdir)) + ".out"
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    with open(outfile, 'w') as f:
        f.write(header)
        while True:
            events = tool.nextEventSet()
            if events is False:
                break
            if events is None:
                continue
            if issubclass(type(events), Event):
                f.write(str(events))
            else:
                for event in events:
                    if issubclass(type(event), list):
                        for ev in event:
                            f.write(str(ev))
                    else:
                        f.write(str(event))


def main():
    parser = argparse.ArgumentParser(description="Project to create unified version of output from different AS tools and"
                                                 "compare them to a ground truth annotation file")

    create_parser = argparse.ArgumentParser(add_help=False)
    create_parser.add_argument("-m", "--majiq_dir", help="directory with 2 majiq files: 1) X.psi.tsv in psi folder 2) output-file of voila tsv run (named: X.voila.tsv)")
    create_parser.add_argument("-s", "--spladder_dir", help="directory of spladder output confirmed.txt files")
    create_parser.add_argument("-w", "--whippet_file", help="whippet-out.psi file")
    create_parser.add_argument("-a", "--asgal_file", help="asgal ASGAL.csv file")
    create_parser.add_argument("-i", "--irfinder_file", help="IRFinder-IR-nondir.txt file")
    create_parser.add_argument("--aspli_dir", help="Directory with three ASpli output files: as_discovery.tab, exon.counts.tab and intron.counts.tab")
    create_parser.add_argument("-e", "--eventpointer_file", help="EventPointer file EventsFound_RNASeq.txt")
    create_parser.add_argument("--sgseq_denovo", help="Tab separated SGSeq file found by de novo analysis, formatted with columns 'from', 'to' and the third one containing the "
                                                      "gene ID and event type.")
    create_parser.add_argument("--sgseq_anno", help="Tab separated SGSeq file based on existing annotation, formatted with columns 'from', 'to' and the third one containing the"
                                                    "gene ID and event type.")

    create_parser.add_argument("-out", "--outdir", help="output directory", required=True)
    create_parser.add_argument("-gtf", "--gtf", help="reference file in gtf format", required=True)
    create_parser.add_argument("-comb", "--combine_me", help="Set this to true if you want MES and MEE events to be counted as ES events "
                                                             "(each skipped exon is one separate ES event)", default=False, action='store_true')
    # create_args = create_parser.parse_args()

    compare_parser = argparse.ArgumentParser(add_help=False)
    compare_parser.add_argument("-a", "--event_annotation", help="Event annotation file for ground truth of events", required=True)
    compare_parser.add_argument("-c", "--compare_file", help="Output of AS tool that will be checked", required=True)
    compare_parser.add_argument("-gtf", "--gtf", help="reference file in gtf format", required=True)
    compare_parser.add_argument("-stats", "--stats_outfile", help="path to outputfile, where statistics will be written to; std-out if not given", required=False)
    compare_parser.add_argument("-comb", "--combine_me", help="Set this to true if you want MES and MEE events to be counted as ES events "
                                                              "(each skipped exon is one separate ES event); should also been used when creating"
                                                              "the compare file!", default=False, action='store_true')
    compare_parser.add_argument("-s", "--strict", help="Use this flag if you want strict comparison between the output and event annotation."
                                                       "Strict means that both start and end coordinate have to be equal so that an event is "
                                                       "counted as correct. Otherwise only one of them has to be equal.", default=False, action='store_true')
    compare_parser.add_argument("-t", "--threshold", help="set threshold to allow for events with minimum distance < threshold to still be"
                                                          "counted as correct; default is 0", default=0)
    # compare_args = compare_parser.parse_args()

    subparser = parser.add_subparsers(help='')
    subparser.add_parser('create', parents=[create_parser], help="Create unified output file(s) for specified input files of AS tools.").set_defaults(func=run_create)
    subparser.add_parser('compare', parents=[compare_parser], help="Compare two files of AS events").set_defaults(func=run_compare)

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
