import gzip
from unified.Event import *

# Type    Interpretation
# CE  Core exon, which may be bounded by one or more alternative AA/AD nodes
# AA  Alternative Acceptor splice site
# AD  Alternative Donor splice site
# RI  Retained intron
# TS  Tandem transcription start site
# TE  Tandem alternative polyadenylation site
# AF  Alternative First exon
# AL  Alternative Last exon
# BS  Circular back-splicing (output only when --circ flag is given)

# init event type dict
from tool_parser.ASParser import ToolParser

EVENT_TYPES = {"CE": "es", "RI": "ir", "AD": "a5", "AA": "a3", "AF": "afe",
               "AL": "ale"}


class WhippetParser(ToolParser):

    NAME = "Whippet"
    current_gene_events = {}                # store all non-NA events of one gene
    current_gene_long_junctions = set()     # store all junctions with more than 1 node inside
    current_gene = None

    def __init__(self, psi_filepath, combine_me):
        super().__init__("Whippet")
        self.combine_me = combine_me                # true -> mes, mee are counted as es events
        self.known_event_types = ["es", "ir", "a3", "a5", "afe", "ale"]
        self.fh = openFile(psi_filepath)
        self.EVENT_PARSING_FUNCTION = {"es": handleES, "ir": handleIR,
                                       "a3": handleA3, "a5": handleA5, "afe": handleAFE,
                                       "ale": handleALE}
        self.header = self.fh.readline().strip().split("\t")
        # dict to save hash of each event -> a way to check for duplicates
        self.event_hash_dict = {EsEvent:[], IrEvent:[], A3Event:[], A5Event:[], MesEvent:[], MeeEvent:[] ,AleEvent:[], AfeEvent:[]}

    def parseLineToEvent(self, event_line):
        event_line = event_line.split("\t")
        gene = event_line[0]
        symbol = event_line[2].split(":")[0]
        if self.current_gene is None: self.current_gene = gene
        node_id = event_line[1]
        idx = gene + ":" + node_id
        strand = event_line[3]
        event_list = None
        if gene != self.current_gene:
            event_list = self.handlePrevGene(strand, self.current_gene, symbol)
            self.current_gene = gene
            self.current_gene_events = {}
            self.current_gene_long_junctions = set()

        # skip NA event types
        if event_line[4] == "NA":
            return event_list
        # skip events with NA for all values
        if sum([x != "NA" for x in event_line[5:14]]) == 0:
            return event_list
        # find out event type
        if event_line[4] in EVENT_TYPES:
            event_type = EVENT_TYPES.get(event_line[4])
        # skip events with unknown event type
        else:
            return event_list

        # coordinate of event
        gen_posi = event_line[2].split(":")
        gen_posi = gen_posi[1].split("-")
        gen_posi = list(map(int, gen_posi))

        # ES event must have exc_edges
        exc_edges = event_line[12]
        if event_type == "es" and exc_edges == "NA":
            return event_list

        # store all junctions with size > 1 to check for MES events later
        all_edges = set([i.split(":")[0] for i in event_line[13].split(",")])
        for edge in all_edges:
            nodes_on_edge = edge.split("-")
            start = int(nodes_on_edge[0])
            end = int(nodes_on_edge[len(nodes_on_edge)-1])   # TODO: i dont know if there are edges in format of: 1-2-3-4; this would be an issue
            if end > start + 2:      # junction covers more than 2 nodes
                self.current_gene_long_junctions.add(edge)

        event = self.handleEvent(gen_posi, event_type, idx, strand, gene, symbol)
        self.current_gene_events[node_id] = event

        return event_list

    # returns None or list of events for one gene
    def nextEventSet(self):
        next_event_line = self.fh.readline()
        if next_event_line == "":
            self.closeFile()
            return False
        next_event_line = next_event_line.strip()
        event_set = self.parseLineToEvent(next_event_line)
        # remove duplicates
        #out_set = set()
        #if event_set is not None:
        #    for event in event_set:
        #        #event_hash = event.__hash__()
        #        if event in self.event_hash_dict[type(event)]:
        #            continue
        #        else:
        #            self.event_hash_dict[type(event)].append(event)
        #            out_set.add(event)
        return event_set

    def handlePrevGene(self, strand, gene, symbol):
        if self.current_gene is None:
            pass
        if len(self.current_gene_long_junctions) == 0:
            return self.current_gene_events.values()
        final_events = set(self.current_gene_events.values())

        # create MES events
        if not self.combine_me:
            for junction in self.current_gene_long_junctions:
                # for each "long" junction check, if more than one ES event is in there -> MES
                j_start = int(junction.split("-")[0])
                j_end = int(junction.split("-")[1])
                es_in_junction = set()              # set of ES events in this junction
                nodes_in_junction = set()           # store node-ids of es events in junction
                for node_id, event in self.current_gene_events.items():
                    # check if this event is inside a "long" junction -> can be part of MES
                    if j_start < int(node_id) < j_end:
                        # check if ES
                        if type(event) is EsEvent:
                            es_in_junction.add(event)
                            nodes_in_junction.add(node_id)

                # more than 1 ES event in this junction -> MES event
                if len(es_in_junction) > 1:
                    idx = self.current_gene + ":" + ":".join(sorted(nodes_in_junction))
                    skipped_exons = [[x.start_skipped, x.end_skipped] for x in es_in_junction]
                    final_events.add(MesEvent(idx, skipped_exons, strand, gene, symbol))
                    final_events.difference_update(es_in_junction)

        return final_events

    # noinspection PyArgumentList
    def handleEvent(self, gen_posi, event_type, idx, strand, gene, symbol):
        if event_type == "FP":
            print("event type is not handled yet. " + idx)
            return None
        func = self.EVENT_PARSING_FUNCTION.get(event_type)
        out = func(gen_posi=gen_posi, strand=str(strand), idx=idx, gene=gene, symbol=symbol)
        return out

    def closeFile(self):
        return self.fh.close()


# exon_skip
def handleES(gen_posi, strand: str, idx, gene, symbol):
    return EsEvent(idx, gen_posi[0], gen_posi[1], strand, gene, symbol)


# intron retention
def handleIR(gen_posi, strand:str, idx, gene, symbol):
    return IrEvent(idx, gen_posi[0], gen_posi[1], strand, gene, symbol)


# alt_3prime
def handleA3(gen_posi, strand:str, idx, gene, symbol):
    #if strand == "+":
    #    gen_posi[1] += 1
    #else:
    #    gen_posi[0] -= 1
    return A3Event(idx, gen_posi[0], gen_posi[1], strand, gene, symbol)


# alt_5prime
def handleA5(gen_posi, strand:str, idx, gene, symbol):
    #if strand == "+":
    #    gen_posi[1] -= 1
    #else:
    #    gen_posi[0] += 1
    return A5Event(idx, gen_posi[0], gen_posi[1], strand, gene, symbol)


def handleAFE(gen_posi, strand:str, idx, gene, symbol):
    return AfeEvent(idx, gen_posi[0], gen_posi[1], strand, gene, symbol)
    #if gen_posi.__eq__(anno_ev.exon1_interval) or gen_posi.__eq__(anno_ev.exon2_interval):
    #    return True
    #return False


def handleALE(gen_posi, strand:str, idx, gene, symbol):
    return AleEvent(idx, gen_posi[0], gen_posi[1], strand, gene, symbol)
    #if gen_posi.__eq__(anno_ev.exon1_interval) or gen_posi.__eq__(anno_ev.exon2_interval):
    #    return True
    #return False


def openFile(filepath):
    if filepath.endswith(".gz"):
        return gzip.open(filepath, mode='rt')
    return open(filepath, mode='r')
