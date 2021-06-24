from gtf_utils.GTFParser import GTFParser
from unified.Event import *


# es: skipped exon start and end
# mes: skipped exons starts and ends
# a3: --: new end -1, actual end; ++: actual start, new start -1
# a5: --: actual start, new start - 1; ++: new end + 1, actual end (index adaption based on gtf)
# mee: 2 lines(mee-->template\ntemplate-->mee): excluded exon start, end
# ir:
# afe:
# ale:

class EventAnnotationReader:

    def __init__(self, event_annotation_file, gtf: GTFParser, combine_me):
        self.combine_me = combine_me                # true -> mes, mee are counted as es events
        self.gtf = gtf
        self.events = {}
        self.count_event_types = {"ES": 0, "MES": 0, "MEE": 0, "IR": 0, "A3": 0, "A5": 0, "AFE": 0, "ALE": 0}
        with open(event_annotation_file, mode='r') as anno:
            header = anno.readline().strip().split("\t")
            idx_number = 0

            for event_line in anno:
                if not event_line.strip():
                    print("Whitespace line; break")
                    break

                event_line = event_line.rstrip('\r\n').split("\t")
                event_type = event_line[0]
                gene = event_line[2].split("_")[0]
                idx = gene + "_" + str(idx_number)
                gtf_gene = gtf.genes.get(gene)
                if gtf_gene is None:
                    continue
                strand = "+" if gtf_gene.is_on_positive_strand() else "-"
                symbol = gtf_gene.chr.id
                start = event_line[3]
                end = event_line[4]

                if gene not in self.events.keys():
                    self.events[gene] = set()

                if event_type == "es":
                    event = {EsEvent(idx, start, end, strand, gene, symbol)}

                elif event_type == "mes":
                    skipped_exons_starts = event_line[3].split(",")
                    if len(skipped_exons_starts) <= 1:
                        continue
                    skipped_exons_ends = event_line[4].split(",")
                    number_skipped_exons = len(skipped_exons_starts)
                    if not self.combine_me:
                        skipped_exons_lists = [[skipped_exons_starts[i], skipped_exons_ends[i]] for i in range(number_skipped_exons)]
                        event = {MesEvent(idx, skipped_exons_lists, strand, gene, symbol)}
                    else:
                        skipped_exons = set([EsEvent(idx, skipped_exons_starts[i], skipped_exons_ends[i], strand, gene, symbol, count=i+1) for i in range(number_skipped_exons)])
                        event = skipped_exons
                        event_type = "es"

                elif event_type == "mee":
                    event2 = anno.readline().strip().split("\t")
                    if not self.combine_me:
                        mee_exons_one = [start, end]
                        mee_exons_two = [event2[3], event2[4]]
                        event = {MeeEvent(idx, [mee_exons_one, mee_exons_two], strand, gene, symbol)}
                    else:
                        event = {EsEvent(idx, start, end, strand, gene, symbol, count=1), EsEvent(idx, event2[3], event2[4], strand, gene, symbol, count=2)}
                        event_type = "es"

                elif event_type == "ir":
                    event = {IrEvent(idx, start, end, strand, gene, symbol)}

                elif event_type == "a3":
                    if strand == "+":
                        event = {A3Event(idx, start, end, strand, gene, symbol)}
                    else:
                        event = {A3Event(idx, start, end, strand, gene, symbol)}

                elif event_type == "a5":
                    if strand == "+":
                        event = {A5Event(idx, start, end, strand, gene, symbol)}
                    else:
                        event = {A5Event(idx, start, end, strand, gene, symbol)}

                elif event_type == "afe":
                    #event2 = anno.readline().strip().split("\t")
                    event = {AfeEvent(idx, start, end, strand, gene, symbol)}

                elif event_type == "ale":
                    #event2 = anno.readline().strip().split("\t")
                    event = {AleEvent(idx, start, end, strand, gene, symbol)}

                else:
                    print("No correct event-type." + idx)
                    continue

                self.events[gene].update(event)
                self.count_event_types[event_type.upper()] += len(event)


                idx_number += 1

