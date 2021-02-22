from event_utils.EventClasses import *
from gtf_utils.GTFParser import GTFParser
from gtf_utils.GeneClass import Gene
from gtf_utils.Interval import Interval
import re
from tool_parser.ASParser import ToolParser

# init event type dict

EVENT_TYPES = {"Cassette Exon": "es", "Mutually Exclusive Exons": "mee", "Retained Intron": "ir",
               "Alternative First Exon": "afe", "Alternative Last Exon": "ale", "Alternative 5' Splice Site": "a5",
               "Alternative 3' Splice Site": "a3", "Complex Event": "FP"}


# exon_skip
def checkES(gene, gen_posi):
    return gen_posi.pseudo_equal(gene.alt_junction)


# intron retention
def checkIR(gene, gen_posi):
    return gene.junction.pseudo_equal(gen_posi)


# mutex_exons
def checkMEE(gene, gen_posi):
    min_exon = min(gene.involved_exons)
    max_exon = max(gene.involved_exons)
    if gene.is_on_positive_strand:
        return gen_posi.pseudo_equal(Interval(gene.exons[min_exon - 1].end, gene.exons[max_exon + 1].start))
    else:
        return gen_posi.pseudo_equal(Interval(gene.exons[max_exon + 1].end, gene.exons[min_exon - 1].start))


# alt_3prime
def checkA3(gene, gen_posi):
    return gen_posi.pseudo_equal(gene.alt_junction)


# alt_5prime
def checkA5(gene, gen_posi):
    return gen_posi.pseudo_equal(gene.alt_junction)


def checkAFE(gene, gen_posi):
    return gen_posi.pseudo_equal(gene.junction) or gen_posi.pseudo_equal(gene.alt_junction)


def checkALE(gene, gen_posi):
    return gen_posi.pseudo_equal(gene.junction) or gen_posi.pseudo_equal(gene.alt_junction)


EVENT_PARSING_FUNCTION = {"es": checkES, "mee": checkMEE, "ir": checkIR,
                          "a3": checkA3, "a5": checkA5, "afe": checkAFE, "ale": checkALE}


class EventPointerParser(ToolParser):

    NAME = "EventPointer"

    def __init__(self, filepath, gtf: GTFParser):
        super().__init__("EventPointer")
        self.fh = self.openFile(filepath)
        self.gtf = gtf
        self.known_event_types = ["es", "mee", "ale", "afe", "ir", "a3", "a5", "FP"]
        self.header = self.fh.readline().strip().split("\t")
        self.header_index = {}
        for i in range(len(self.header)):
            self.header_index[self.header[i]] = i

    def parseLineToEvent(self, event_line):
        event_line = event_line.split("\t")
        event_type = event_line[self.header_index["Event Type"]]
        if event_type in EVENT_TYPES:
            event_type = EVENT_TYPES.get(event_type)
        else:
            event_type = "FP"
        gene = event_line[self.header_index["Gene"]]
        # skip if gene is unknwon
        if gene in self.gtf.genes:
            gene = self.gtf.genes.get(gene)
        else:
            return None, None, None, True
        # skip if event type of gene is not known
        if len(gene.event_types) > 0 and gene.event_types[0] not in self.known_event_types:
            return None, None, None, True
        gen_posi = event_line[self.header_index["Genomic Position"]].split(":")
        chr = gen_posi[0]
        gen_posi = gen_posi[1].split("-")
        gen_posi = Interval(gen_posi[0], gen_posi[1])
        return self.checkEvent(gene, gen_posi, event_type)

    def checkEvent(self, gene: Gene, gen_posi: Interval, event_type: str):
        if len(gene.event_types) == 0:
            return gene, event_type, False, False
        if event_type == "FP":
            return gene, event_type, False, False
        if event_type != gene.event_types[0]:
            return gene, event_type, False, False
        func = EVENT_PARSING_FUNCTION.get(event_type)
        ev_correct = func(gene, gen_posi)
        return gene, event_type, ev_correct, False

    def nextEvent(self):
        next_event_line = self.fh.readline()
        if next_event_line == "":
            self.closeFile()
            return None, None, None, False
        next_event_line = next_event_line.strip()
        return self.parseLineToEvent(next_event_line)

    def openFile(self, filepath):
        return open(filepath, mode='r')

    def closeFile(self):
        return self.fh.close()
