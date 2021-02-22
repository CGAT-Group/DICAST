import re

from gtf_utils.GTFParser import GTFParser
from gtf_utils.Interval import Interval
from tool_parser.ASParser import ToolParser

# init event type dict
EVENT_TYPES = {"SE": "es", "RI": "ir", "MXE": "mee", "A5SS": "a5", "A3SS": "a3", "AFE": "afe",
               "ALE": "ale"}


# exon_skip
def checkES(events):
    gene, transcript, is_on_positive_strand, gen_posi = events[0][1:5]
    return gene.alt_junction.pseudo_equal(gen_posi) or gene.alt_counter_junction.pseudo_equal(gen_posi)


# intron retention
def checkIR(events):
    gene, transcript, is_on_positive_strand, gen_posi = events[0][1:5]
    return gene.junction.pseudo_equal(gen_posi)


# mutex_exons
def checkMEE(events):
    # chromosome, gene, transcripts, exon_pre_interval, exon1_interval, exon2_interval, exon_aft_interval
    gene, transcript, is_on_positive_strand, gen_posi = events[0][1:5]
    min_exon = min(gene.involved_exons)
    max_exon = max(gene.involved_exons)
    if is_on_positive_strand:
        return gen_posi.pseudo_equal(Interval(gene.exons[min_exon - 1].end, gene.exons[max_exon + 1].start))
    else:
        return gen_posi.pseudo_equal(Interval(gene.exons[max_exon + 1].end, gene.exons[min_exon - 1].start))


# alt_3prime
def checkA3(events):
    # chromosome, gene, transcripts, event_start, event_end
    gene, transcript, is_on_positive_strand, gen_posi = events[0][1:5]
    return gen_posi.pseudo_equal(gene.alt_junction)


# alt_5prime
def checkA5(events):
    # chromosome, gene, transcripts, event_start, event_end
    gene, transcript, is_on_positive_strand, gen_posi = events[0][1:5]
    return gen_posi.pseudo_equal(gene.alt_junction)


def checkAFE(events):
    gene, transcript, is_on_positive_strand, gen_posi = events[0][1:5]
    transcript2 = events[1][2]
    gen_posi2 = events[1][4]
    skipped_exon1 = gene.exons[gene.involved_exons[0]]
    skipped_exon2 = gene.exons[gene.involved_exons[1]]
    if is_on_positive_strand:
        ends_correct = Interval(skipped_exon1.start, skipped_exon2.start).pseudo_equal(Interval(gen_posi.start, gen_posi2.start))
        exon_start_correct = gene.exons[max(gene.involved_exons) + 1].start == gen_posi.end
        return ends_correct and exon_start_correct
    ends_correct = Interval(skipped_exon1.end, skipped_exon2.end).pseudo_equal(Interval(gen_posi.start, gen_posi2.start))
    exon_end_correct = gene.exons[max(gene.involved_exons) + 1].end == gen_posi.end
    return ends_correct and exon_end_correct


def checkALE(events):
    gene, transcript, is_on_positive_strand, gen_posi = events[0][1:5]
    transcript2 = events[1][2]
    gen_posi2 = events[1][4]
    skipped_exon1 = gene.exons[gene.involved_exons[0]]
    skipped_exon2 = gene.exons[gene.involved_exons[1]]
    if is_on_positive_strand:
        ends_correct = Interval(skipped_exon1.end, skipped_exon2.end).pseudo_equal(Interval(gen_posi.end, gen_posi2.end))
        exon_start_correct = gene.exons[min(gene.involved_exons) - 1].end == gen_posi.start
        return ends_correct and exon_start_correct
    ends_correct = Interval(skipped_exon1.start, skipped_exon2.start).pseudo_equal(Interval(gen_posi.end, gen_posi2.end))
    exon_end_correct = gene.exons[min(gene.involved_exons) - 1].start == gen_posi.start
    return ends_correct and exon_end_correct


EVENT_PARSING_FUNCTION = {"es": checkES, "mee": checkMEE, "ir": checkIR,
                          "a3": checkA3, "a5": checkA5, "afe": checkAFE, "ale": checkALE}


class SGSeqDenovoReader(ToolParser):

    NAME = "SGSeq"

    def __init__(self, filepath, gtf: GTFParser):
        super().__init__("SGSeq_denovo")
        self.gtf = gtf
        self.fh = self.openFile(filepath)
        self.known_event_types = ["es", "mee", "ir", "a3", "a5", "afe", "ale", "FP"]
        self.header = self.fh.readline().rstrip().replace("\"", "").split(",")
        self.header_index = {}
        for i in range(len(self.header)):
            self.header_index[self.header[i]] = i
        self.reg = re.compile(",(?!\\s)(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)")

    # collects corresponding events (1/2, 2/2, ...)
    def collectEvents(self, event_line):
        chr, gene, transcript, is_on_positive_strand, gen_posi, event_type, current_event_index, corresponding_events_number, val = self.parseLineToEvent(
            event_line)
        events = [[] for i in range(corresponding_events_number)]
        events[0] = [chr, gene, transcript, is_on_positive_strand, gen_posi, event_type, current_event_index,
                     corresponding_events_number, val]
        for i in range(1, corresponding_events_number):
            next_event_line = self.fh.readline().rstrip('\r\n')
            if next_event_line.endswith(", "):
                next_event_line += self.fh.readline().rstrip('\r\n')
            events[i] = list(self.parseLineToEvent(next_event_line))
        if gene not in self.gtf.genes or event_line[-1] == "NA":
            return None, None, None, True
        gene = self.gtf.genes.get(gene)
        for arr in events:
            arr[1] = gene
        if len(gene.event_types) > 0 and gene.event_types[0] not in self.known_event_types:
            return None, None, None, True
        if len(gene.event_types) == 0:
            return gene, event_type, False, False
        if event_type == "FP":
            return gene, event_type, False, False
        if corresponding_events_number > 2:
            return gene, "FP", False, False
        if gene.event_types[0] != event_type:
            return gene, event_type, False, False
        #chr = gene.chr
        #transcript = self.gtf.transcripts[transcript]
        func = EVENT_PARSING_FUNCTION.get(event_type)
        ev_correct = func(events)
        return gene, event_type, ev_correct, False

    def parseLineToEvent(self, event_line):
        event_line = self.reg.split(event_line)
        #print(event_line)
        event_type = self.parseEventType(event_line[self.header_index["variantType"]])
        # collect events belonging to each other:
        event_index = event_line[self.header_index["variantName"]].split("/")
        current_event_index = int(event_index[0][-1])
        corresponding_events_number = int(event_index[1][0])

        start = event_line[self.header_index["from"]].split(":")
        end = event_line[self.header_index["to"]].split(":")
        gen_posi = Interval(start[2], end[2])
        chr = start[1]
        gene = event_line[self.header_index["geneName"]]
        transcript = event_line[self.header_index["txName"]]
        is_on_positive_strand = start[3] == "+"
        val = event_line[self.header_index["test"]]
        return chr, gene, transcript, is_on_positive_strand, gen_posi, event_type, current_event_index, corresponding_events_number, val

    def nextEvent(self):
        next_event_line = self.fh.readline()
        if next_event_line == "":
            self.closeFile()
            return None, None, None, False
        next_event_line = next_event_line.rstrip('\r\n')
        if next_event_line.endswith(", "):
            next_event_line += self.fh.readline().rstrip('\r\n')
        return self.collectEvents(next_event_line)

    def openFile(self, filepath):
        return open(filepath, mode='r')

    def closeFile(self):
        return self.fh.close()

    def parseEventType(self, variant_type):
        if variant_type.startswith("c("):
            return "FP"
        variant_type = variant_type.split(":")[0]
        if variant_type not in EVENT_TYPES:
            return "FP"
        return EVENT_TYPES[variant_type]
