from unified.Event import IrEvent

# init event type dict
from tool_parser.ASParser import ToolParser


def openFile(filepath):
    return open(filepath, mode='r')


EVENT_TYPES = {"IR": "ir"}


class IRFinderParser(ToolParser):

    NAME = "IRFinder"

    def __init__(self, filepath):
        super().__init__("IRFinder")
        self.fh = openFile(filepath)
        self.header = self.fh.readline().strip().split("\t")
        self.known_event_types = ["ir"]

    def parseLineToEvent(self, event_line):
        event_line = event_line.split("\t")
        start = int(event_line[1]) + 1
        end = int(event_line[2])
        gene = event_line[3].split("/")[0]
        ir_ratio = float(event_line[19])
        if ir_ratio == 0:
            return None
        strand = event_line[5]
        symbol = event_line[0]
        event = IrEvent(idx=gene, intron_start=start, intron_end=end, strand=strand, gene=gene, symbol=symbol)
        return event

    def nextEventSet(self):
        next_event_line = self.fh.readline()
        if next_event_line == "":
            self.closeFile()
            return False
        next_event_line = next_event_line.strip()
        return self.parseLineToEvent(next_event_line)

    def closeFile(self):
        return self.fh.close()



