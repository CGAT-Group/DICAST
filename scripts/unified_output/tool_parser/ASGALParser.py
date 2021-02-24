from unified.Event import *
from gtf_utils.GeneClass import Gene
from gtf_utils.GTFParser import GTFParser
from tool_parser.ASParser import ToolParser


def openFile(filepath):
    return open(filepath, mode='r')


# init event type dict
EVENT_TYPES = {"ES": "es", "IR": "ir", "A3": "a3", "A5": "a5"}


class ASGALParser(ToolParser):

    NAME = "ASGAL"

    def __init__(self, filepath, gtf: GTFParser, combine_me):
        super().__init__("ASGAL")
        self.combine_me = combine_me            # true -> mes, mee are counted as es events
        self.gtf = gtf
        self.fh = openFile(filepath)
        self.header = self.fh.readline().strip().split("\t")
        self.known_event_types = ["es", "mes", "ir", "a3", "a5"]
        self.EVENT_PARSING_FUNCTION = {"es": checkES, "ir": checkIR,
                                       "a3": checkA3A5, "a5": checkA3A5}

    def nextEventSet(self):
        next_event_line = self.fh.readline()
        if next_event_line == "":
            self.closeFile()
            return False
        next_event_line = next_event_line.strip()
        return self.parseLineToEvent(next_event_line)

    def parseLineToEvent(self, event_line):
        event_line = event_line.split(",")
        event_type = EVENT_TYPES.get(event_line[0])
        if event_type is None:
            return None
        junction = (event_line[1], event_line[2])
        transcripts = event_line[4].split("/")
        gene = transcripts[0].split("_")[0]
        gene = self.gtf.genes.get(gene)
        chromosome = gene.chr.id
        idx = str(gene.feature_id + "_" + str(event_line[1]) + "_" + str(event_line[2]))

        if gene.is_on_positive_strand():
            strand = "+"
        else:
            strand = "-"
        event = self.checkEvent(junction, gene, event_type, strand, idx, chromosome)

        return event

    def checkEvent(self, junction, gene: Gene, event_type, strand, idx, symbol):
        func = self.EVENT_PARSING_FUNCTION.get(event_type)
        return func(event_type, junction, gene, strand, idx, symbol, combine_me=self.combine_me)

    def closeFile(self):
        return self.fh.close()


# exon_skip
def checkES(event_type, junction, gene: Gene, strand, idx, symbol, combine_me):
    j_start = int(junction[0])
    j_end = int(junction[1])
    exons_in_gtf = gene.getExonsSortedByPosition()
    skipped_exons = []

    # for each exon for this gene, check if it is in the ES junction
    for exon in exons_in_gtf:
        if exon.end < j_start or exon.start > j_end:
            continue
        # check if current exon is in junction
        if exon.start >= j_start and exon.end <= j_end:
            skipped_exons.append([str(exon.start), str(exon.end)])

    if len(skipped_exons) == 0:
        print("Found no exons that fit into this ES junction. " + idx)
        return None
    if combine_me:
        skipped_exons = set([EsEvent(idx, skipped_exons[i][0], skipped_exons[i][1], strand, gene.feature_id, symbol, count=i+1) for i in
                             range(len(skipped_exons))])
        return skipped_exons
    else:
        if len(skipped_exons) == 1:
            return {EsEvent(idx, skipped_exons[0][0], skipped_exons[0][1], strand, gene.feature_id, symbol)}
        else:
            return {MesEvent(idx, skipped_exons, strand, gene.feature_id, symbol)}


# intron retention
def checkIR(event_type, junction, gene: Gene, strand, idx, symbol, combine_me):
    return {IrEvent(idx, junction[0], junction[1], strand, gene.feature_id, symbol)}


# alt_3 and alt_5 prime
# get alternative part of event by looking at exon and junction (intron) boundaries
# NOTE: this does not handle events, where a junction covers more than 1 exon!
def checkA3A5(event_type, junction, gene: Gene, strand, idx, symbol, combine_me):
    j_start = int(junction[0])
    j_end = int(junction[1])
    exons_in_gtf = gene.getExonsSortedByPosition()
    event = None
    count=1

    # iterate over all exons for this gene (until second to last)
    for e in range(len(exons_in_gtf)-1):
        exon = exons_in_gtf[e]
        next_exon = exons_in_gtf[e+1]

        if (event_type == "a3" and strand == "-") or (event_type == "a5" and strand == "+"):
            if j_end+1 == next_exon.start:
                alt_part = (str(exon.end), str(j_start))
                if j_start < exon.end:                          # check if junction is not inside of exon
                    alt_part = (str(j_start-1), str(exon.end))
                    count=999
                #if not (exon.end <= j_start <= next_exon.start):  # check if junction really starts in between two exons; it might cover more than 1 exon
                #    break
                if j_start == exon.end:                         # alt_part start & stop can be the same -> NaN for alt_start
                    alt_part = ("nan", str(j_start))
                if event_type == "a3":
                    return {A3Event(idx, alt_part[0], alt_part[1], strand, gene.feature_id, symbol, count)}
                else:
                    return {A5Event(idx, alt_part[0], alt_part[1], strand, gene.feature_id, symbol, count)}

        elif (event_type == "a5" and strand == "-") or (event_type == "a3" and strand == "+"):
            if j_start-1 == exon.end:
                alt_part = (str(j_end), str(next_exon.start))
                if next_exon.start < j_end:                             # check if junction is not inside of exon
                    alt_part = (str(next_exon.start), str(j_end+1))
                    count=999
                #if not (exon.end <= j_end <= next_exon.start):            # check if junction really ends in between two exons; it might cover more than 1 exon
                #    break
                if j_end == next_exon.start:                            # alt_part start & stop can be the same -> NaN for alt_start
                    alt_part = (str(j_end), "nan")
                if event_type == "a3":
                    return {A3Event(idx, alt_part[0], alt_part[1], strand, gene.feature_id, symbol, count)}
                else:
                    return {A5Event(idx, alt_part[0], alt_part[1], strand, gene.feature_id, symbol, count)}



    return {event}

