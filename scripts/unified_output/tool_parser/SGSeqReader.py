from gtf_utils.GTFParser import GTFParser
from gtf_utils.Interval import Interval
from tool_parser.ASParser import ToolParser
from itertools import product, chain
from . import EventPointerParser
from unified.Event import IrEvent, EsEvent, MeeEvent
from gtf_utils.ExonClass import Exon

# init event type dict
EVENT_TYPES = {"SE": "es", "RI": "ir", "MXE": "mee", "A5SS": "a5", "A3SS": "a3", "AFE": "afe",
               "ALE": "ale"}


# exon_skip
def checkES(events):
    result = set()
    for event in events:
        chrm, gene, strand, gen_posi, event_type, current_event_index, corresponding_events_number, combine_me = event
        if gene is None or event_type != "es":
            continue
        if current_event_index == 1:
            if strand == "+":
                result.add(
                    EventPointerParser.createES(gene=gene.feature_id, chrm=chrm, start=gen_posi.start, end=gen_posi.end, strand=strand, gtf_gene=gene, combine_me=combine_me))
            else:
                skipped = None
                for transcript in gene.transcripts.values():
                    exons = transcript.getExons()
                    for index, exon in enumerate(exons):
                        if not exon.getEnd() == gen_posi.end:
                            continue
                        try:
                            if exons[index + 2].getStart() == gen_posi.start:
                                skipped = exons[index + 1]
                                break
                        except IndexError:
                            continue

                    if skipped is not None:
                        break

                if skipped is None:
                    skipped = Exon(None, None, None, "nan", "nan", None, None, None)

                es_event = EsEvent(idx=gene.feature_id,
                                   start_skipped=skipped.getStart(),
                                   end_skipped=skipped.getEnd(),
                                   strand=strand,
                                   gene=gene.feature_id,
                                   symbol=chrm)
                result.add(es_event)

    return result


# intron retention
def checkIR(events):
    result = set()
    for event in events:
        chrm, gene, strand, gen_posi, event_type, current_event_index, corresponding_events_number, combine_me = event
        if gene is None or event_type != "ir":
            continue
        if current_event_index == 1:
            if strand == "+":
                start = gen_posi.start + 1
                end = gen_posi.end - 1
            else:
                start = gen_posi.end + 1
                end = gen_posi.start - 1

            ir_event = IrEvent(idx=gene.feature_id,
                               intron_start=start,
                               intron_end=end,
                               strand=strand,
                               gene=gene.feature_id,
                               symbol=chrm)

            result.add(ir_event)

    return result


# mutex_exons
def checkMEE(events):
    result = []
    for event in events:
        chrm, gene, strand, gen_posi, event_type, current_event_index, corresponding_events_number, combine_me = event
        if gene is None or event_type != "mee":
            continue
        if current_event_index == 1:
            # if exon.getEnd() != end or exons[index + 2].getStart() != start:

            if strand == "+":
                result.append(
                    EventPointerParser.createMEE(gene=gene.feature_id, chrm=chrm, start=gen_posi.start, end=gen_posi.end, strand=strand, gtf_gene=gene, combine_me=combine_me))
            else:
                exon_1_start = "nan"
                exon_2_start = "nan"
                exon_1_end = "nan"
                exon_2_end = "nan"
                found_alt_1 = False
                found_alt_2 = False
                for transcript in gene.transcripts.values():
                    exons = transcript.getExons()
                    for index, exon in enumerate(exons):
                        try:
                            if exon.getEnd() != gen_posi.end or exons[index + 2].getStart() != gen_posi.start:
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
                    es_events = [EsEvent(gene.feature_id, exon_1_start, exon_1_end, strand, gene.feature_id, chrm, count=1),
                                 EsEvent(gene.feature_id, exon_2_start, exon_2_end, strand, gene.feature_id, chrm, count=2)]
                    result.append(es_events)
                else:
                    mee_event = MeeEvent(idx=gene.feature_id,
                                         mee_exons=[(exon_1_start, exon_1_end),
                                                    (exon_2_start, exon_2_end)],
                                         strand=strand,
                                         gene=gene.feature_id,
                                         symbol=chrm)
                    result.append(mee_event)

    return result


# alt_3prime
def checkA3(events):
    result = set()
    for event in events:
        chrm, gene, strand, gen_posi, event_type, current_event_index, corresponding_events_number, combine_me = event
        if gene is None or event_type != "a3":
            continue
        if current_event_index == 1:
            result.add(EventPointerParser.createA3(gene=gene.feature_id, chrm=chrm, start=gen_posi.start, end=gen_posi.end, strand=strand, gtf_gene=gene, combine_me=combine_me))
    return result


# alt_5prime
def checkA5(events):
    result = set()
    for event in events:
        chrm, gene, strand, gen_posi, event_type, current_event_index, corresponding_events_number, combine_me = event
        if gene is None or event_type != "a5":
            continue
        if current_event_index == 1:
            result.add(EventPointerParser.createA5(gene=gene.feature_id, chrm=chrm, start=gen_posi.start, end=gen_posi.end, strand=strand, gtf_gene=gene, combine_me=combine_me))
    return result


def checkAFE(events):
    result = set()
    for event in events:
        chrm, gene, strand, gen_posi, event_type, current_event_index, corresponding_events_number, combine_me = event
        if gene is None or event_type != "afe" or gen_posi is None:
            continue
        result.add(EventPointerParser.createAFE(gene=gene.feature_id, chrm=chrm, start=gen_posi.start, end=gen_posi.end, strand=strand, gtf_gene=gene, combine_me=combine_me))

    return result


def checkALE(events):
    result = set()
    for event in events:
        chrm, gene, strand, gen_posi, event_type, current_event_index, corresponding_events_number, combine_me = event
        if gene is None or event_type != "ale" or gen_posi is None:
            continue
        result.add(EventPointerParser.createALE(gene=gene.feature_id, chrm=chrm, start=gen_posi.start, end=gen_posi.end, strand=strand, gtf_gene=gene, combine_me=combine_me))

    return result


EVENT_PARSING_FUNCTION = {"es": checkES, "mee": checkMEE, "ir": checkIR,
                          "a3": checkA3, "a5": checkA5, "afe": checkAFE, "ale": checkALE}


class SGSeqReader(ToolParser):
    NAME = "SGSeq_denovo"

    def __init__(self, filepath, gtf: GTFParser, combine_me):
        super().__init__("SGSeq")
        self.gtf = gtf
        self.fh = self.openFile(filepath)
        self.known_event_types = ["es", "mee", "ir", "a3", "a5", "afe", "ale", "FP"]
        self.header = self.fh.readline().rstrip().split("\t")
        self.header_index = {}
        for i in range(len(self.header)):
            self.header_index[self.header[i].strip()] = i
        self.combine_me = combine_me

    # collects corresponding events (1/2, 2/2, ...)
    def collectEvents(self, event_line):
        chrom, genes, is_on_positive_strand, gen_posi, event_types, current_event_index, corresponding_events_number = self.parseLineToEvent(
            event_line)
        if genes is None and event_types is None:
            return None

        all_possible_events = list(product(genes, event_types))
        events = []
        for event_combination in all_possible_events:
            events.append([chrom, self.gtf.genes.get(event_combination[0]), is_on_positive_strand, gen_posi, event_combination[1], current_event_index,
                           corresponding_events_number, self.combine_me])

        for i in range(0, corresponding_events_number - 1):
            next_event_line = self.fh.readline().rstrip('\r\n')
            _, _, _, new_gen_posi, _, new_event_index, _ = self.parseLineToEvent(next_event_line)
            for event_combination in all_possible_events:
                events.append([chrom, self.gtf.genes.get(event_combination[0]), is_on_positive_strand, new_gen_posi, event_combination[1], new_event_index,
                               corresponding_events_number, self.combine_me])

        _events = [EVENT_PARSING_FUNCTION.get(event)(events) for event in EVENT_TYPES.values()]
        events_flatten = list(chain(*_events))
        return events_flatten

    def parseLineToEvent(self, full_event_line):
        # event_line = self.reg.split(full_event_line)
        event_line = full_event_line.split("\t")

        event_types = self.parseEventTypes(event_line[2])
        if len(event_types) == 0:
            return None, None, None, None, None, None, None

        # collect events belonging to each other:
        event_indices = event_line[2].split("/")
        current_event_index = int(event_indices[0][-1])
        corresponding_events_number = int(event_indices[1][0])

        start_column_split = event_line[self.header_index["from"]].split(":")
        end_column_split = event_line[self.header_index["to"]].split(":")
        gen_posi = Interval(start_column_split[2].strip(), end_column_split[2].strip())
        chr = start_column_split[1].strip()
        genes = event_line[2].split("_")[0].split(",")
        genes = [gene.strip() for gene in genes]
        strand = start_column_split[3].strip()

        return chr, genes, strand, gen_posi, event_types, current_event_index, corresponding_events_number

    def nextEventSet(self):
        next_event_line = self.fh.readline()
        if next_event_line == "":
            self.closeFile()
            return False
        next_event_line = next_event_line.rstrip('\r\n')
        if next_event_line.endswith(", "):
            next_event_line += self.fh.readline().rstrip('\r\n')
        return self.collectEvents(next_event_line)

    def openFile(self, filepath):
        return open(filepath, mode='r')

    def closeFile(self):
        return self.fh.close()

    def parseEventTypes(self, type_column):
        event_types = type_column.strip().split("_")[-1].split(",")
        recognized_event_types = []
        for event_type in event_types:
            if event_type in EVENT_TYPES:
                recognized_event_types.append(event_type)
        return [EVENT_TYPES[event_type] for event_type in recognized_event_types]
