from unified.EventAnnotationReader import EventAnnotationReader
from event_utils.EventClasses import *
from gtf_utils.GeneClass import Gene
from gtf_utils.Interval import Interval

# init event type dict
from tool_parser.ASParser import ToolParser


def openFile(filepath):
    return open(filepath, mode='r')


EVENT_TYPES = {"IR": "ir"}


class IRFinder(ToolParser):

    def __init__(self, filepath, gtf, event_anno: EventAnnotationReader):
        super().__init__("IRFinder")
        self.gtf = gtf
        self.event_anno = event_anno
        self.fh = openFile(filepath)
        self.header = self.fh.readline().strip().split("\t")
        self.known_event_types = ["ir"]
        self.EVENT_PARSING_FUNCTION = {"ir": self.checkIR}

    def parseLineToEvent(self, event_line):
        event_line = event_line.split("\t")
        start = int(event_line[1]) + 1
        end = int(event_line[2])
        skipped_intron = Interval(start, end)
        gene = event_line[3].split("/")[0]
        event_type = "ir"
        if gene not in self.gtf.genes:
            return None, None, None, True
        gene = self.gtf.genes.get(gene)
        if len(gene.event_types) == 1 and gene.event_types[0] not in self.known_event_types:
            return None, None, None, True
        ev_correct = self.checkEvent(gene, skipped_intron, event_type)
        return gene, event_type, ev_correct, False

    def nextEvent(self):
        next_event_line = self.fh.readline()
        if next_event_line == "":
            self.closeFile()
            return None, None, None, None
        next_event_line = next_event_line.strip()
        return self.parseLineToEvent(next_event_line)

    def closeFile(self):
        return self.fh.close()

    def checkEvent(self, gene: Gene, gen_posi: Interval, event_type):
        if len(gene.event_types) == 0:
            return False
        anno_ev = self.event_anno.events.get(gene.getID())
        func = self.EVENT_PARSING_FUNCTION.get(event_type)
        return func(gene, gen_posi, anno_ev)

    # intron retention
    def checkIR(self, gene, gen_posi: Interval, anno_ev: IREvent):
        if gen_posi.pseudo_equal(anno_ev.intron_interval):
            return True
        return False
