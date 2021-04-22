from unified.Event import IrEvent, A3Event, A5Event, EsEvent, AfeEvent, AleEvent, MeeEvent
from gtf_utils.GTFParser import GTFParser
from gtf_utils.GeneClass import Gene
from gtf_utils.ExonClass import Exon
from gtf_utils.Interval import Interval
import re
from tool_parser.ASParser import ToolParser

# init event type dict

EVENT_TYPES = {"Cassette Exon": "es", "Mutually Exclusive Exons": "mee", "Retained Intron": "ir",
               "Alternative First Exon": "afe", "Alternative Last Exon": "ale", "Alternative 5' Splice Site": "a5",
               "Alternative 3' Splice Site": "a3"}


# exon_skip
def createES(gene, chrm, start, end, strand, gtf_gene: Gene, combine_me):
    skipped = None

    for transcript in gtf_gene.transcripts.values():
        exons = transcript.getExons()
        for index, exon in enumerate(exons):
            if not exon.getEnd() == start:
                continue
            try:
                if exons[index + 2].getStart() == end:
                    skipped = exons[index + 1]
                    break
            except IndexError:
                continue

        if skipped is not None:
            break

    if skipped is None:
        skipped = Exon(None, None, None, "nan", "nan", None, None, None)

    return EsEvent(idx=gene,
                   start_skipped=skipped.getStart(),
                   end_skipped=skipped.getEnd(),
                   strand=strand,
                   gene=gene,
                   symbol=chrm)


# intron retention
def createIR(gene, chrm, start, end, strand, gtf_gene: Gene, combine_me):
    return IrEvent(idx=gene,
                   intron_start=start + 1,
                   intron_end=end - 1,
                   strand=strand,
                   gene=gene,
                   symbol=chrm)


# mutex_exons
def createMEE(gene, chrm, start, end, strand, gtf_gene: Gene, combine_me):
    exon_1_start = "nan"
    exon_2_start = "nan"
    exon_1_end = "nan"
    exon_2_end = "nan"
    found_alt_1 = False
    found_alt_2 = False
    for transcript in gtf_gene.transcripts.values():
        exons = transcript.getExons()
        for index, exon in enumerate(exons):
            try:
                if exon.getEnd() != start or exons[index + 2].getStart() != end:
                    continue
            except IndexError:
                continue

            middle_exon_start = exons[index + 1].getStart()
            middle_exon_end = exons[index + 1].getEnd()

            if not found_alt_1:
                exon_1_start = middle_exon_start
                exon_1_end = middle_exon_end
                found_alt_1 = True

            if not found_alt_2 and middle_exon_start != exon_1_start and middle_exon_end != exon_1_end:
                exon_2_start = middle_exon_start
                exon_2_end = middle_exon_end
                found_alt_2 = True


        if found_alt_1 and found_alt_2:
            break

    if combine_me:
        return [EsEvent(gene, exon_1_start, exon_1_end, strand, gene, chrm, count=1),
                EsEvent(gene, exon_2_start, exon_2_end, strand, gene, chrm, count=2)]
    else:
        return MeeEvent(idx=gene,
                        mee_exons=[(exon_1_start, exon_1_end),
                                   (exon_2_start, exon_2_end)],
                        strand=strand,
                        gene=gene,
                        symbol=chrm)


# alt_3prime
def createA3(gene, chrm, start, end, strand, gtf_gene: Gene, combine_me):
    exon_start = "nan"
    exon_end = "nan"
    found_start_exon = False
    found_end_exon = False
    for transcript in gtf_gene.transcripts.values():
        exons = transcript.getExons()
        for index, exon in enumerate(exons):
            if strand == "+":
                if exon.getEnd() == start:
                    if not found_start_exon:
                        try:
                            exon_start = exons[index + 1].getStart()
                        except IndexError:
                            continue
                        found_start_exon = True

                if exon.getStart() == end:
                    exon_end = end - 1
                    found_end_exon = True
                    break
            else:
                if exon.getStart() == start:
                    if exons[index - 1].getEnd() == end:
                        continue
                    try:
                        exon_end = exons[index - 1].getEnd()
                    except IndexError:
                        continue
                    found_start_exon = True
                if exon.getEnd() == end:
                    exon_start = end + 1
                    found_end_exon = True
                    break

        if found_start_exon and found_end_exon:
            break

    return A3Event(idx=gene,
                   alt_start=exon_start,
                   alt_end=exon_end,
                   strand=strand,
                   gene=gene,
                   symbol=chrm)


# alt_5prime
def createA5(gene, chrm, start, end, strand, gtf_gene: Gene, combine_me):
    exon_start = "nan"
    exon_end = "nan"
    found_start_exon = False
    found_end_exon = False
    for transcript in gtf_gene.transcripts.values():
        exons = transcript.getExons()
        for index, exon in enumerate(exons):
            if strand == "+":
                if exon.getStart() == end:
                    if not found_start_exon:
                        try:
                            exon_end = exons[index - 1].getEnd()
                        except IndexError:
                            continue
                        found_start_exon = True
                if exon.getEnd() == start:
                    exon_start = start + 1
                    found_end_exon = True
                    break
            else:
                if exon.getEnd() == end:
                    if exons[index + 1].getStart() == start:
                        continue
                    try:
                        exon_start = exons[index + 1].getStart()
                    except IndexError:
                        continue
                    found_start_exon = True
                if exon.getStart() == start:
                    exon_end = start - 1
                    found_end_exon = True
                    break

        if found_start_exon and found_end_exon:
            break

    return A5Event(idx=gene,
                   alt_start=exon_start,
                   alt_end=exon_end,
                   strand=strand,
                   gene=gene,
                   symbol=chrm)


# TODO: AFE and ALE don't work, the reported coordinates are nonsense
def createAFE(gene, chrm, start, end, strand, gtf_gene: Gene, combine_me):
    return AfeEvent(idx=gene,
                    afe_start=start,
                    afe_end=end,
                    strand=strand,
                    gene=gene,
                    symbol=chrm)


def createALE(gene, chrm, start, end, strand, gtf_gene: Gene, combine_me):
    return AleEvent(idx=gene,
                    ale_start=start,
                    ale_end=end,
                    strand=strand,
                    gene=gene,
                    symbol=chrm)


EVENT_PARSING_FUNCTION = {"es": createES, "mee": createMEE, "ir": createIR,
                          "a3": createA3, "a5": createA5, "afe": createAFE, "ale": createALE}


class EventPointerParser(ToolParser):

    NAME = "EventPointer"

    def __init__(self, filepath, gtf: GTFParser, combine_me):
        super().__init__("EventPointer")
        self.fh = self.openFile(filepath)
        self.gtf = gtf
        self.known_event_types = ["es", "mee", "ale", "afe", "ir", "a3", "a5"]
        self.header = self.fh.readline().strip().split("\t")
        self.header_index = {}
        for i in range(len(self.header)):
            self.header_index[self.header[i]] = i
        self.combine_me = combine_me

    def parseLineToEvent(self, event_line):
        event_line = event_line.split("\t")
        event_type = event_line[self.header_index["Event Type"]]
        if event_type in EVENT_TYPES:
            event_type = EVENT_TYPES.get(event_type)
        else:
            return None

        gene = event_line[self.header_index["Gene"]]
        # skip if gene is unknown
        if gene in self.gtf.genes:
            gtf_gene = self.gtf.genes.get(gene)
        else:
            return None

        if gtf_gene.is_on_positive_strand():
            strand = "+"
        else:
            strand = "-"

        gene_positions = event_line[self.header_index["Genomic Position"]].split(":")
        chrm = gene_positions[0]
        gene_positions = gene_positions[1].split("-")
        start = int(gene_positions[0])
        end = int(gene_positions[1])

        func = EVENT_PARSING_FUNCTION.get(event_type)
        event = func(gene, chrm, start, end, strand, gtf_gene, self.combine_me)
        return event

    def nextEventSet(self):
        next_event_line = self.fh.readline()
        if next_event_line == "":
            self.closeFile()
            return False
        next_event_line = next_event_line.strip()
        return self.parseLineToEvent(next_event_line)

    def openFile(self, filepath):
        return open(filepath, mode='r')

    def closeFile(self):
        return self.fh.close()
