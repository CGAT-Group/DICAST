from event_utils.EventClasses import *
from gtf_utils.ExonClass import Exon
from gtf_utils.GTFParser import GTFParser
from gtf_utils.GeneClass import Gene
from gtf_utils.Interval import Interval
import re
from tool_parser.ASParser import ToolParser

# init event type dict

EVENT_TYPES = {"ES": "es", "IR": "ir", "Alt5ss": "a5", "Alt3ss": "a3"}


# exon_skip
def checkES(gene: Gene, exon_number, junc1, junc2, junc3):
    junc1 = parseJunction(junc1)
    junc2 = parseJunction(junc2)
    junc3 = parseJunction(junc3)
    #    if junc1 is None or junc2 is None or junc3 is None:
    if junc1 is None and junc2 is None and junc3 is None:
        return False

    if gene.is_on_positive_strand():
        junc = Interval(gene.exons.get(exon_number).end, gene.exons.get(exon_number + 1).start)
    else:
        junc = Interval(gene.exons.get(exon_number + 1).end, gene.exons.get(exon_number).start)
    if junc1 is not None:
        if not (junc1.pseudo_equal(gene.junction) or junc1.pseudo_equal(junc)):
            return False
    if junc2 is not None:
        if not (junc2.pseudo_equal(junc) or junc2.pseudo_equal(gene.junction)):
            return False
    if junc3 is not None:
        if not (junc3.pseudo_equal(gene.alt_junction) or junc3.pseudo_equal(gene.alt_junction)):
            return False
    else:
        return False
    return True


#    if gene.is_on_positive_strand():
#        junc = Interval(gene.exons.get(exon_number).end, gene.exons.get(exon_number + 1).start)
#        return junc1.pseudo_equal(gene.junction) and junc2.pseudo_equal(junc) and junc3.pseudo_equal(gene.alt_junction)
#    else:
#        junc = Interval(gene.exons.get(exon_number + 1).end, gene.exons.get(exon_number).start)
#        return junc1.pseudo_equal(junc) and junc2.pseudo_equal(gene.junction) and junc3.pseudo_equal(gene.alt_junction)


# intron retention
def checkIR(gene: Gene, exon_number, junc1, junc2, junc3):
    junc = parseJunction(junc3)
    if junc is None:
        return None
    return junc.pseudo_equal(gene.junction)


# alt_3prime
def checkA3(gene: Gene, exon_number, junc1, junc2, junc3):
    junc1 = parseJunction(junc1)
    junc2 = parseJunction(junc2)
    junc3 = parseJunction(junc3)
    if junc3 is None:
        return False
    if not junc3.pseudo_equal(gene.alt_junction):
        return False
    if junc1 is not None:
        if not gene.junction.pseudo_equal(junc1):
            return False
    if junc2 is not None:
        if not gene.junction.pseudo_equal(junc2):
            return False
    # if sum([junc1 is None, junc2 is None, junc3 is None]) > 1:
    #     return False
    # if junc1 is None:
    #     junc = junc2
    #     alt_junc = junc3
    # elif junc2 is None:
    #     junc = junc1
    #     alt_junc = junc3
    # elif junc3 is None:
    #     junc = junc1
    #     alt_junc = junc2
    # else:
    #     print("A3 has 3 junctions!!!")
    #     exit(1)
    # ev_correct = gene.junction.pseudo_equal(junc) and gene.alt_junction.pseudo_equal(alt_junc)
    # ev_correct = ev_correct or (gene.junction.pseudo_equal(alt_junc) and gene.alt_junction.pseudo_equal(junc))
    # return ev_correct
    return True


# alt_5prime
def checkA5(gene: Gene, exon_number, junc1, junc2, junc3):
    junc1 = parseJunction(junc1)
    junc2 = parseJunction(junc2)
    junc3 = parseJunction(junc3)
    if junc3 is None:
        return False
    if not junc3.pseudo_equal(gene.alt_junction):
        return False
    if junc1 is not None:
        if not gene.junction.pseudo_equal(junc1):
            return False
    if junc2 is not None:
        if not gene.junction.pseudo_equal(junc2):
            return False
    # if sum([junc1 is None, junc2 is None, junc3 is None]) > 1:
    #     return False
    # if junc1 is None:
    #     junc = junc2
    #     alt_junc = junc3
    # elif junc2 is None:
    #     junc = junc1
    #     alt_junc = junc3
    # elif junc3 is None:
    #     junc = junc1
    #     alt_junc = junc2
    # else:
    #     print("A5 has 3 junctions!!!")
    #     exit(1)
    # ev_correct = gene.junction.pseudo_equal(junc) and gene.alt_junction.pseudo_equal(alt_junc)
    # ev_correct = ev_correct or (gene.junction.pseudo_equal(alt_junc) and gene.alt_junction.pseudo_equal(junc))
    # return ev_correct
    return True


def parseJunction(junction_string: str):
    if junction_string == "NA":
        return None
    if junction_string.__contains__(";"):
        return None
    chr, start, stop = junction_string.split(".")
    return Interval(start, stop)


EVENT_PARSING_FUNCTION = {"es": checkES, "ir": checkIR, "a3": checkA3, "a5": checkA5}


class ASpliParser(ToolParser):

    NAME = "ASpli"

    def __init__(self, filepath, gtf: GTFParser):
        super().__init__("ASpli")
        self.gtf = gtf
        self.known_event_types = ["es", "ir", "a3", "a5", "FP"]
        self.fh = self.openFile(filepath)
        self.header = self.fh.readline().rstrip().split("\t")
        self.header_index = {}
        for i in range(len(self.header)):
            self.header_index[self.header[i]] = i

    def parseLineToEvent(self, event_line):
        event_line = event_line.split("\t")
        event_type = event_line[self.header_index["event"]]
        # skip NA events (--> no evidence, only parsed from gtf)
        if sum([x != 'NA' for x in event_line[2:]]) == 0:
            return None, None, None, True
        if event_type == "-":
            return None, None, None, True
        if event_type in EVENT_TYPES:
            event_type = EVENT_TYPES.get(event_type)
        else:
            event_type = "FP"
        gene, exon = event_line[self.header_index[""]].split(":")
        exon = exon.split("_")[0]
        if gene in self.gtf.genes:
            gene = self.gtf.genes.get(gene)
        else:
            return None, None, None, True
        if len(gene.event_types) > 0 and gene.event_types[0] not in self.known_event_types:
            return None, None, None, True
        # parse exon number from Aspli to gtf format
        exon_number: int = self.getExonFromAspli(gene, exon)
        if exon_number is None:
            return None, None, None, True
        # parse junctions
        junc1 = event_line[self.header_index["J1"]]
        junc2 = event_line[self.header_index["J2"]]
        junc3 = event_line[self.header_index["J3"]]
        return self.checkEvent(gene, exon_number, event_type, junc1, junc2, junc3)

    def checkEvent(self, gene: Gene, exon_number: int, event_type: str, junc1: str, junc2: str, junc3: str):
        if len(gene.event_types) == 0:
            return gene, event_type, False, False
        if event_type == "FP":
            return gene, event_type, False, False
        if gene.event_types[0] != event_type:
            return gene, event_type, False, False
        func = EVENT_PARSING_FUNCTION.get(event_type)
        ev_correct = func(gene, exon_number, junc1, junc2, junc3)
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

    def getExonFromAspli(self, gene: Gene, exon: str):
        # delete leading e and zeros
        exon_number = exon[1:].lstrip("0")
        try:
            exon_number = int(exon_number)
        except (ValueError):
            return None
        # handle negative strand (annotation in aspli is reversed: E001 == highest exon number in gtf)
        if not gene.is_on_positive_strand():
            exon_number = len(gene.exons) - exon_number + 1
        return exon_number
