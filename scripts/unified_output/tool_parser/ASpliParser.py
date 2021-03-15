from gtf_utils.ExonClass import Exon
from gtf_utils.GTFParser import GTFParser
from gtf_utils.GeneClass import Gene
from gtf_utils.Interval import Interval
import re
from tool_parser.ASParser import ToolParser
from unified.Event import IrEvent, A3Event, A5Event, EsEvent
import pandas as pd
import sys

# init event type dict
EVENT_TYPES = {"ES": "es", "IR": "ir", "Alt5ss": "a5", "Alt3ss": "a3"}


# exon_skip
def createES(gene: Gene, start, end, strand):
    return EsEvent(idx=gene.feature_id,
                   start_skipped=start,
                   end_skipped=end,
                   strand=strand,
                   gene=gene.feature_id,
                   symbol=gene.chr.id)


# intron retention
def createIR(gene: Gene, start, end, strand):
    return IrEvent(idx=gene.feature_id,
                   intron_start=start,
                   intron_end=end,
                   strand=strand,
                   gene=gene.feature_id,
                   symbol=gene.chr.id)


# alt_3prime
def createA3(gene: Gene, start, end, strand):
    if strand == "-":
        start = start - 1
        end = end - 1

    return A3Event(idx=gene.feature_id,
                   alt_start=start,
                   alt_end=end+1,
                   strand=strand,
                   gene=gene.feature_id,
                   symbol=gene.chr.id)


# alt_5prime
def createA5(gene: Gene, start, end, strand):
    if strand == "+":
        start = start - 1
        end = end - 1
    return A5Event(idx=gene.feature_id,
                   alt_start=start,
                   alt_end=end+1,
                   strand=strand,
                   gene=gene.feature_id,
                   symbol=gene.chr.id)


EVENT_PARSING_FUNCTION = {"es": createES, "ir": createIR, "a3": createA3, "a5": createA5}


class ASpliParser(ToolParser):

    NAME = "ASpli"

    def __init__(self, discovery_file, exon_counts_file, intron_counts_file, gtf: GTFParser):
        super().__init__("ASpli")
        self.gtf = gtf
        self.discovery_file_handler = self.openFile(discovery_file)
        self.discovery_file_header = self.discovery_file_handler.readline().rstrip().split("\t")
        self.discovery_file_header_index = {}
        for i in range(len(self.discovery_file_header)):
            self.discovery_file_header_index[self.discovery_file_header[i]] = i

        # Get exon and intron counts files as data frames
        # First columns, containing the "Gene:Exon" event ID, have no column name
        self.exon_counts_df = pd.read_csv(exon_counts_file, header=0, usecols=["Unnamed: 0", "start", "end"], delimiter="\t")
        self.exon_counts_df.rename(columns={"Unnamed: 0": "event_id"}, inplace=True)
        self.exon_counts_df.set_index("event_id", inplace=True)

        self.intron_counts_df = pd.read_csv(intron_counts_file, header=0, usecols=["Unnamed: 0", "start", "end"], delimiter="\t")
        self.intron_counts_df.rename(columns={"Unnamed: 0": "event_id"}, inplace=True)
        self.intron_counts_df.set_index("event_id", inplace=True)

    def parseLineToEvent(self, event_line):
        event_line = event_line.split("\t")
        event_type = event_line[self.discovery_file_header_index["event"]]
        bam_columns = event_line[-2:]

        # Handle * events same as normal events
        event_type = event_type.replace("*", "")

        # Skip event if a BAM column is NA, or the event type is not recognized
        if 'NA' in bam_columns or event_type not in EVENT_TYPES:
            return None

        event_type = EVENT_TYPES.get(event_type)
        event_id = event_line[self.discovery_file_header_index[""]]
        gene, _ = event_id.split(":")
        if gene in self.gtf.genes:
            gene = self.gtf.genes.get(gene)
        else:
            print("Gene not found GTF file. Wrong file?")
            return None

        if gene.is_on_positive_strand():
            strand = "+"
        else:
            strand = "-"

        if event_type == "ir":
            row = self.intron_counts_df.loc[event_id]
        else:
            row = self.exon_counts_df.loc[event_id]

        func = EVENT_PARSING_FUNCTION.get(event_type)
        event = func(gene, row["start"], row["end"], strand)
        return event

    def nextEventSet(self):
        next_event_line = self.discovery_file_handler.readline()
        if next_event_line == "":
            self.closeFile()
            return False
        next_event_line = next_event_line.strip()
        return self.parseLineToEvent(next_event_line)

    def openFile(self, filepath):
        try:
            return open(filepath, mode='r')
        except FileNotFoundError:
            print("Couldn't find file %s, please make sure that the path is correct and the file is properly named." % filepath)
            sys.exit(1)

    def closeFile(self):
        self.discovery_file_handler.close()

