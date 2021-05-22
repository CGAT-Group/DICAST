from gtf_utils.GenomicRegionClass import GenomicRegion
from gtf_utils.Interval import Interval


class Gene(GenomicRegion):
    def __init__(self, feature_id, chr, start, end, is_on_positive_strand):
        super().__init__(feature_id, start, end, is_on_positive_strand)
        self.chr = chr
        self.exons = {}
        self.transcripts = {}
        self.event_types = list()
        # counter: how often seen in event prediction
        self.seen = 0
        self.correct_seen = False
        self.junction: Interval = None
        self.counter_junction: Interval = None
        self.alt_junction: Interval = None
        self.alt_counter_junction: Interval = None
        self.involved_exons = None

    def getChromosome(self):
        return self.chr

    def addTranscript(self, transcript):
        self.transcripts[transcript.getID()] = transcript
        tr_event_types = transcript.getEventTypes()
        if not len(tr_event_types) == 0:
            for ev_type in tr_event_types:
                self.event_types.append(ev_type)

    def addExon(self, exon):
        self.exons[exon.getID()] = exon

    def getExonsSortedByPosition(self):
        return [val for key, val in sorted(self.exons.items(), reverse=(not self.is_on_positive_strand()))]

    def getExonsSortedByNumber(self):
        return [val for key, val in sorted(self.exons.items())]

    # !!! watch out !!!
    # es & mes & ir working on dicts and gene numbers as key
    # rest working with indices
    # counter is for MAJIQ: source and target handling
    def calcEventJunctions(self):
        if len(self.transcripts) == 1:
            return
        for tr in self.transcripts.values():
            if tr.is_template:
                template_tr = tr
            else:
                event_tr = tr
        ev_type = event_tr.event_types[0]
        if ev_type == "es":
            skipped_ex_number = set([ex.getID() for ex in template_tr.exons]).symmetric_difference(
                set([ex.getID() for ex in event_tr.exons])).pop()
            junc_counter = None
            alt_junc_counter = None
            if self.is_on_positive_strand():
                junc = Interval(self.exons[skipped_ex_number - 1].end, self.exons[skipped_ex_number].start)
                alt_junc = Interval(self.exons[skipped_ex_number - 1].end, self.exons[skipped_ex_number + 1].start)
                if len(self.exons) > skipped_ex_number > 1:
                    junc_counter = Interval(self.exons[skipped_ex_number - 1].end, self.exons[skipped_ex_number].start)
                    alt_junc_counter = Interval(self.exons[skipped_ex_number].end,
                                                self.exons[skipped_ex_number + 1].start)
            else:
                junc = Interval(self.exons[skipped_ex_number].end, self.exons[skipped_ex_number - 1].start)
                alt_junc = Interval(self.exons[skipped_ex_number + 1].end, self.exons[skipped_ex_number - 1].start)
                if len(self.exons) > skipped_ex_number > 1:
                    junc_counter = Interval(self.exons[skipped_ex_number + 1].end, self.exons[skipped_ex_number].start)
                    alt_junc_counter = Interval(self.exons[skipped_ex_number + 1].end,
                                                self.exons[skipped_ex_number - 1].start)
            self.junction, self.alt_junction = junc, alt_junc
            self.counter_junction, self.alt_counter_junction = junc_counter, alt_junc_counter
            self.involved_exons = [skipped_ex_number]
        elif ev_type == "mes":
            skipped_ex_numbers = set([ex.getID() for ex in template_tr.exons]).symmetric_difference(
                set([ex.getID() for ex in event_tr.exons]))
            skipped_ex_numbers = list(skipped_ex_numbers)
            min_skipped_ex_number = min(skipped_ex_numbers)
            max_skipped_ex_number = max(skipped_ex_numbers)
            junc_counter = None
            alt_junc_counter = None
            if self.is_on_positive_strand():
                junc = Interval(self.exons[min_skipped_ex_number - 1].end, self.exons[min_skipped_ex_number].start)
                alt_junc = Interval(self.exons[min_skipped_ex_number - 1].end,
                                    self.exons[max_skipped_ex_number + 1].start)
                if len(self.exons) > max_skipped_ex_number and min_skipped_ex_number > 1:
                    junc_counter = Interval(self.exons[min_skipped_ex_number - 1].end,
                                            self.exons[max_skipped_ex_number + 1].start)
                    alt_junc_counter = Interval(self.exons[max_skipped_ex_number].end,
                                                self.exons[max_skipped_ex_number + 1].start)
            else:
                junc = Interval(self.exons[min_skipped_ex_number].end, self.exons[min_skipped_ex_number - 1].start)
                alt_junc = Interval(self.exons[max_skipped_ex_number + 1].end,
                                    self.exons[min_skipped_ex_number - 1].start)
                if len(self.exons) > max_skipped_ex_number and min_skipped_ex_number > 1:
                    junc_counter = Interval(self.exons[max_skipped_ex_number + 1].end,
                                            self.exons[min_skipped_ex_number - 1].start)
                    alt_junc_counter = Interval(self.exons[max_skipped_ex_number + 1].end,
                                                self.exons[max_skipped_ex_number].start)
            self.junction, self.alt_junction = junc, alt_junc
            self.counter_junction, self.alt_counter_junction = junc_counter, alt_junc_counter
            self.involved_exons = [skipped_ex_numbers]
        elif ev_type == "mee":
            skipped_ex_numbers = set([ex.getID() for ex in template_tr.exons]).symmetric_difference(
                set([ex.getID() for ex in event_tr.exons]))
            self.involved_exons = list(skipped_ex_numbers)
            return
            # self.junction, self.alt_junction = junc, alt_junc
        elif ev_type == "a3":
            x = -1
            if self.is_on_positive_strand():
                for i in range(len(self.exons)):
                    if template_tr.exons[i].start != event_tr.exons[i].start:
                        x = i
                        break
                junc = Interval(template_tr.exons[x - 1].end, template_tr.exons[x].start)
                alt_junc = Interval(event_tr.exons[x - 1].end, event_tr.exons[x].start)
            else:
                for i in range(len(self.exons)):
                    if template_tr.exons[i].end != event_tr.exons[i].end:
                        x = i
                        break
                junc = Interval(template_tr.exons[x].end, template_tr.exons[x + 1].start)
                alt_junc = Interval(event_tr.exons[x].end, event_tr.exons[x + 1].start)
            self.junction, self.alt_junction = junc, alt_junc
            self.involved_exons = [template_tr.exons[x].feature_id]
        elif ev_type == "a5":
            x = -1
            if self.is_on_positive_strand():
                for i in range(len(self.exons)):
                    if template_tr.exons[i].end != event_tr.exons[i].end:
                        x = i
                        break
                junc = Interval(template_tr.exons[x].end, template_tr.exons[x + 1].start)
                alt_junc = Interval(event_tr.exons[x].end, event_tr.exons[x + 1].start)
            else:
                for i in range(len(self.exons)):
                    if template_tr.exons[i].start != event_tr.exons[i].start:
                        x = i
                        break
                junc = Interval(template_tr.exons[x - 1].end, template_tr.exons[x].start)
                alt_junc = Interval(event_tr.exons[x - 1].end, event_tr.exons[x].start)
            self.junction, self.alt_junction = junc, alt_junc
            self.involved_exons = [template_tr.exons[x].feature_id]
        elif ev_type == "afe":
            skipped_ex_numbers = set([ex.getID() for ex in template_tr.exons]).symmetric_difference(
                set([ex.getID() for ex in event_tr.exons]))
            self.involved_exons = list(skipped_ex_numbers)
            if self.is_on_positive_strand():
                junc = Interval(self.exons[max(skipped_ex_numbers)].end, self.exons[max(skipped_ex_numbers) + 1].start)
                alt_junc = Interval(self.exons[min(skipped_ex_numbers)].end,
                                    self.exons[max(skipped_ex_numbers) + 1].start)
            else:
                junc = Interval(self.exons[max(skipped_ex_numbers) + 1].end, self.exons[max(skipped_ex_numbers)].start)
                alt_junc = Interval(self.exons[max(skipped_ex_numbers) + 1].end,
                                    self.exons[min(skipped_ex_numbers)].start)
            self.junction, self.alt_junction = junc, alt_junc
        elif ev_type == "ale":
            skipped_ex_numbers = set([ex.getID() for ex in template_tr.exons]).symmetric_difference(
                set([ex.getID() for ex in event_tr.exons]))
            self.involved_exons = list(skipped_ex_numbers)
            if self.is_on_positive_strand():
                junc = Interval(self.exons[min(skipped_ex_numbers) - 1].end, self.exons[min(skipped_ex_numbers)].start)
                alt_junc = Interval(self.exons[min(skipped_ex_numbers) - 1].end,
                                    self.exons[max(skipped_ex_numbers)].start)
            else:
                junc = Interval(self.exons[max(skipped_ex_numbers)].end, self.exons[min(skipped_ex_numbers) - 1].start)
                alt_junc = Interval(self.exons[min(skipped_ex_numbers)].end,
                                    self.exons[min(skipped_ex_numbers) - 1].start)
            self.junction, self.alt_junction = junc, alt_junc
        elif ev_type == "ir":
            missing_exon_number = set([ex.getID() for ex in template_tr.exons]).symmetric_difference(
                set([ex.getID() for ex in event_tr.exons])).pop()
            # remember: exon_numbers != index; exons sorted by gen_posi (reverse of exon_number)
            if self.is_on_positive_strand():
                junc = Interval(self.exons[missing_exon_number - 1].end, self.exons[missing_exon_number].start)
            else:
                junc = Interval(self.exons[missing_exon_number].end, self.exons[missing_exon_number - 1].start)
            self.junction, self.alt_junction = junc, None
