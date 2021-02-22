from unified.Event import *
import os
import gzip
from tool_parser.ASParser import ToolParser

EVENT_TYPES = {"exon_skip": "es", "mutex_exons": "mee", "intron_retention": "ir", "mult_exon_skip": "mes",
               "alt_5prime": "a5", "alt_3prime": "a3"}


# exon_skip
def parseES(es_line, combine_me):
    es_line = es_line.split("\t")
    gene = es_line[3]
    symbol = es_line[0]
    idx = gene + es_line[2]
    strand = es_line[1]

    return {EsEvent(idx, es_line[6], es_line[7], strand, gene, symbol)}


# mult_exon_skip
def parseMES(mes_line, combine_me):
    mes_line = mes_line.split("\t")
    strand = mes_line[1]
    gene = mes_line[3]
    symbol = mes_line[0]
    idx = gene + mes_line[2]
    skipped_exons_starts = mes_line[6].split(":")
    skipped_exons_ends = mes_line[7].split(":")
    number_skipped_exons = len(skipped_exons_starts)
    skipped_exons_tuples = [(skipped_exons_starts[i], skipped_exons_ends[i]) for i in
                               range(number_skipped_exons)]
    if combine_me:
        skipped_exons = set([EsEvent(idx, skipped_exons_tuples[i][0], skipped_exons_tuples[i][1], strand, gene, symbol, count=i+1) for i in
                             range(len(skipped_exons_tuples))])
        return skipped_exons
    else:
        return {MesEvent(idx, skipped_exons_tuples, strand, gene, symbol)}


# intron retention
def parseIR(ir_line, combine_me):
    ir_line = ir_line.split("\t")
    strand = ir_line[1]
    gene = ir_line[3]
    symbol = ir_line[0]
    idx = gene + ir_line[2]

    return {IrEvent(idx, ir_line[6], ir_line[7], strand, gene, symbol)}


# mutex_exons
def parseMEE(mee_line, combine_me):
    mee_line = mee_line.split("\t")
    strand = mee_line[1]
    gene = mee_line[3]
    symbol = mee_line[0]
    idx = gene + mee_line[2]
    mutex1 = (mee_line[6], mee_line[7])
    mutex2 = (mee_line[8], mee_line[9])
    mee_lst = [mutex1, mutex2]

    if combine_me:
        out = {EsEvent(idx, mutex1[0], mutex1[1], strand, gene, symbol, count=1), EsEvent(idx, mutex2[0], mutex2[1], strand, gene, symbol, count=2)}
        return out
    else:
        return {MeeEvent(idx, mee_lst, strand, gene, symbol)}


# alt_3prime
def parseA3(alt_3prime_line, combine_me):
    alt_3prime_line = alt_3prime_line.split("\t")
    gene = alt_3prime_line[3]
    symbol = alt_3prime_line[0]
    idx = gene + alt_3prime_line[2]
    strand = alt_3prime_line[1]

    if strand == "+":
        return {A3Event(idx, alt_3prime_line[8], alt_3prime_line[6], strand, gene, symbol)}
    else:
        return {A3Event(idx, alt_3prime_line[7], alt_3prime_line[9], strand, gene, symbol)}


# alt_5prime
def parseA5(alt_5prime_line, combine_me):
    alt_5prime_line = alt_5prime_line.split("\t")
    gene = alt_5prime_line[3]
    symbol = alt_5prime_line[0]
    idx = gene + alt_5prime_line[2]
    strand = alt_5prime_line[1]

    if strand == "+":
        return {A5Event(idx, alt_5prime_line[7], alt_5prime_line[9], strand, gene, symbol)}
    else:
        return {A5Event(idx, alt_5prime_line[8], alt_5prime_line[6], strand, gene, symbol)}


EVENT_PARSING_FUNCTION = {"es": parseES, "mes": parseMES, "mee": parseMEE, "ir": parseIR,
                          "a3": parseA3, "a5": parseA5}



class SplAdderParser(ToolParser):

    NAME = "SplAdder"

    def __init__(self, spladder_output_directory, combine_me):

        super().__init__("SplAdder")
        self.combine_me = combine_me
        self.input_files = list()
        self.event_types_found = list()
        self.event_type_done = list()
        self.current_pointer = 0
        self.current_event_type = None
        self.fh = None

        for filename in os.listdir(spladder_output_directory):
            if filename.endswith(".confirmed.txt"):
                event_type_found = \
                    filename.replace("merge_graphs_", "").rsplit("_", 1)[
                        0]
                event_type_found = EVENT_TYPES[event_type_found]
                self.event_types_found.append(event_type_found)
                input_file = os.path.join(spladder_output_directory, filename)
                self.input_files.append(input_file)
                print("adding event_type: " + event_type_found + " (" + input_file + ")")


        if len(self.event_types_found) == 0:
            print("NO EVENTS FOUND!!!")
            return
        else:
            self.current_event_type = self.event_types_found[self.current_pointer]
            # init first filehandle
            self.fh = self.openFile(self.input_files[self.current_pointer])
            # skip header
            self.fh.readline()

    def nextEventSet(self):
        event = None
        if self.fh is None:
            return False
        next_event_line = self.fh.readline().strip()
        if next_event_line == "":
            print("no more events of type: " + self.current_event_type)
            self.event_type_done.append(self.current_event_type)
            self.current_pointer += 1
            if self.current_pointer >= len(self.event_types_found):
                print("no more events found. finished")
                self.closeFile()
                return False
            else:
                self.current_event_type = self.event_types_found[self.current_pointer]
                print("continue with events of type: " + self.current_event_type)
                self.closeFile()
                self.fh = self.openFile(self.input_files[self.current_pointer])
                # skip header
                self.fh.readline()
                return self.nextEventSet()
        return self.parseLineToEvent(next_event_line, self.current_event_type)


    def parseLineToEvent(self, event_line, event_type):

        func = EVENT_PARSING_FUNCTION.get(event_type)
        event = func(event_line, self.combine_me)

        return event


    def openFile(self, filepath):
        if filepath.endswith(".gz"):
            return gzip.open(filepath, mode='rt')
        return open(filepath, mode='r')

    def closeFile(self):
        return self.fh.close()
