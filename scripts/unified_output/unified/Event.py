# this function can be used to sort skipped exons by coordinate value (for example for MES events);
# inputs is:
# list of skipped exons; each index of pattern: "start-end"
# one of these coordinates can be "nan" -> these exons are not sorted, but will be added to the end of the list
# return:
# single list with same pattern, sorted by coordinate
def sort_exons_with_na(skipped_exons):
    starts = []
    ends = []
    nan_pairs = []
    for i, exon in enumerate(skipped_exons):
        start = exon[0]
        end = exon[1]
        if start == "nan" or end == "nan":
            nan_pairs.append((start, end))
        else:
            starts.append(start)
            ends.append(end)

    # the remaining coordinates can just be sorted like this; the exon-pairs will keep the same index
    sorted_starts = sorted(starts, key=int)
    sorted_ends = sorted(ends, key=int)

    # add the nan pair to the end of each list
    for pair in nan_pairs:
        sorted_starts.append(pair[0])
        sorted_ends.append(pair[1])

    sorted_list = [[sorted_starts[i], sorted_ends[i]] for i in range(len(sorted_starts))]
    return sorted_list


class Event:

    id = None
    start = 0
    end = 0
    gene = None
    symbol = None
    exon_list = None
    count = 1

    def __init__(self, idx, gene, symbol, count):
        self.id = idx + ":" + str(count)
        #self.id = idx
        self.gene = gene
        self.symbol = symbol
        self.count = count

    def to_string(self):
        print("empty event.")

    def get_start_stop(self):
        return [self.start, self.end]



class EsEvent(Event):

    start_skipped = 0
    end_skipped = 0
    strand = ""
    count = 0

    def __init__(self, idx,  start_skipped, end_skipped, strand, gene, symbol, count=1):
        super().__init__(idx, gene, symbol, count)
        self.start_skipped = start_skipped
        self.end_skipped = end_skipped
        self.strand = strand
        self.count = count

    def __eq__(self, other):
        if type(self) != type(other): return False
        return self.start_skipped == other.start_skipped and self.end_skipped == other.end_skipped

    def __hash__(self):
        return hash((self.gene, self.start_skipped, self.end_skipped))

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.symbol, self.gene, self.id, self.strand, "ES", self.count, self.start_skipped, self.end_skipped)

    def to_string(self):
        print("ES-Event: " + str(self.id) + "  Start skipped: " + str(self.start_skipped) + "; End skipped: " + str(self.end_skipped))

    def get_start_stop(self):
        return [self.start_skipped, self.end_skipped]


class MeeEvent(Event):
    strand = ""
    mee_exons = []
    count = 0

    def __init__(self, idx, mee_exons, strand, gene, symbol, count=1):
        super().__init__(idx, gene, symbol, count)
        self.mee_exons = sort_exons_with_na(mee_exons)
        self.strand = strand
        self.symbol = symbol
        self.count = count

    def to_string(self):
        print("MEE-Event: " + str(self.id) + " mut. exclusive exons are: " + str(self.mee_exons))

    def __str__(self):
        mee_exons_starts = ""
        mee_exons_ends = ""
        for i in range(len(self.mee_exons)-1):
            mee_exons_starts += str(self.mee_exons[i][0]) + ","
            mee_exons_ends += str(self.mee_exons[i][1]) + ","
        mee_exons_starts += str(self.mee_exons[len(self.mee_exons)-1][0])
        mee_exons_ends += str(self.mee_exons[len(self.mee_exons)-1][1])
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.symbol, self.gene, self.id, self.strand, "MEE", self.count, mee_exons_starts, mee_exons_ends)

    def __hash__(self):
        return hash((self.gene, frozenset(frozenset(ev) for ev in self.mee_exons)))

    def __eq__(self, other):
        if type(self) != type(other): return False
        if len(self.mee_exons) != len(other.mee_exons): return False
        b = True
        for i in range(len(self.mee_exons)):
            if self.mee_exons[i] is not other.mee_exons[i]:
                b = False
        return b

    def get_start_stop(self):
        return self.mee_exons


class MesEvent(Event):

    mes_skipped_exons = []
    strand = ""
    count = 0

    def __init__(self, idx,  mes_skipped, strand, gene, symbol, count=1):
        super().__init__(idx, gene, symbol, count)
        self.mes_skipped_exons = sort_exons_with_na(mes_skipped)
        self.strand = strand
        self.symbol = symbol
        self.count = count

    def __eq__(self, other):
        if type(self) != type(other): return False
        if len(self.mes_skipped_exons) != len(other.mes_skipped_exons): return False
        b = True
        for i in range(len(self.mes_skipped_exons)):
            if self.mes_skipped_exons[i] is not other.mes_skipped_exons[i]:
                b = False
        return b

    def __hash__(self):
        return hash((self.gene, frozenset(frozenset(ev) for ev in self.mes_skipped_exons)))

    def __str__(self):
        skipped_exons_starts = ""
        skipped_exons_ends = ""
        for i in range(len(self.mes_skipped_exons)-1):
            skipped_exons_starts += str(self.mes_skipped_exons[i][0]) + ","
            skipped_exons_ends += str(self.mes_skipped_exons[i][1]) + ","
        skipped_exons_starts += str(self.mes_skipped_exons[len(self.mes_skipped_exons)-1][0])
        skipped_exons_ends += str(self.mes_skipped_exons[len(self.mes_skipped_exons)-1][1])
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.symbol, self.gene, self.id, self.strand, "MES", self.count, skipped_exons_starts, skipped_exons_ends)

    def to_string(self):
        print("MES-Event: " + str(self.id) + " Skipped exons are: " + str(self.mes_skipped_exons))

    def get_start_stop(self):
        return self.mes_skipped_exons


class A3Event(Event):

    alt_start = 0
    alt_end = 0
    strand = ""
    count = 0

    def __init__(self, idx, alt_start, alt_end, strand, gene, symbol, count=1):
        super().__init__(idx, gene, symbol, count)
        self.alt_end = alt_end
        self.alt_start = alt_start
        self.strand = strand
        self.count = count
        self.symbol = symbol

    def __eq__(self, other):
        if type(self) != type(other): return False
        return self.alt_start == other.alt_start and self.alt_end == other.alt_end

    def __hash__(self):
        return hash((self.gene, self.alt_start, self.alt_end))

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.symbol, self.gene, self.id, self.strand, "A3", self.count, self.alt_start, self.alt_end)

    def to_string(self):
        print("A3-Event: " + str(self.id) + "  Start alternative: " + str(self.alt_start) + "; End alternative: " + str(self.alt_end))

    def get_start_stop(self):
        return [self.alt_start, self.alt_end]


class A5Event(Event):
    alt_start = 0
    alt_end = 0
    strand = ""
    count = 0

    def __init__(self, idx, alt_start, alt_end, strand, gene, symbol, count=1):
        super().__init__(idx, gene, symbol, count)
        self.alt_end = alt_end
        self.alt_start = alt_start
        self.strand = strand
        self.count = count
        self.symbol = symbol

    def __eq__(self, other):
        if type(self) != type(other): return False
        return self.alt_start == other.alt_start and self.alt_end == other.alt_end

    def __hash__(self):
        return hash((self.gene, self.alt_start, self.alt_end))

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.symbol, self.gene, self.id, self.strand, "A5", self.count, self.alt_start, self.alt_end)

    def to_string(self):
        print("A5-Event: " + str(self.id) + "  Start alternative: " + str(self.alt_start) + "; End alternative: " + str(self.alt_end))

    def get_start_stop(self):
        return [self.alt_start, self.alt_end]


class IrEvent(Event):

    intron_start = 0
    intron_end = 0
    strand = ""
    count = 0

    def __init__(self, idx, intron_start, intron_end, strand, gene, symbol, count=1):
        super().__init__(idx, gene, symbol, count)
        self.intron_start = intron_start
        self.intron_end = intron_end
        self.strand = strand
        self.count = count
        self.symbol = symbol

    def __eq__(self, other):
        if type(self) != type(other): return False
        return self.intron_start == other.intron_start and self.intron_end == other.intron_end

    def __hash__(self):
        return hash((self.gene, self.intron_start, self.intron_end))

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.symbol, self.gene, self.id, self.strand, "IR", self.count, self.intron_start, self.intron_end)

    def to_string(self):
        print("IR-Event: " + str(self.id) + "  Start intron: " + str(self.intron_start) + "; End intron: " + str(self.intron_end))

    def get_start_stop(self):
        return [self.intron_start, self.intron_end]


class AfeEvent(Event):

    afe_start = 0
    afe_end = 0
    strand = ""
    count = 0

    def __init__(self, idx, afe_start, afe_end, strand, gene, symbol, count=1):
        super().__init__(idx, gene, symbol, count)
        self.afe_start = afe_start
        self.afe_end = afe_end
        self.strand = strand
        self.symbol = symbol
        self.count = count

    def to_string(self):
        print("AFE-Event: " + str(self.id) + " Start AFE: " + str(self.afe_start) + "; End AFE: " + str(self.afe_end))

    def __hash__(self):
        return hash((self.gene, self.afe_start, self.afe_end))

    def __eq__(self, other):
        if type(self) != type(other): return False
        return self.afe_start == other.afe_start and self.afe_end == other.afe_end

    def get_start_stop(self):
        return[self.afe_start, self.afe_end]

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.symbol, self.gene, self.id, self.strand, "AFE", self.count,
                                                     self.afe_start, self.afe_end)


class AleEvent(Event):
    ale_start = 0
    ale_end = 0
    strand = ""
    count = 0

    def __init__(self, idx, ale_start, ale_end, strand, gene, symbol, count=1):
        super().__init__(idx, gene, symbol, count)
        self.ale_start = ale_start
        self.ale_end = ale_end
        self.strand = strand
        self.count = count

    def to_string(self):
        print("ALE-Event: " + str(self.id) + " Start AFE: " + str(self.ale_start) + "; End AFE: " + str(self.ale_end))

    def __hash__(self):
        return hash((self.gene, self.ale_start, self.ale_end))

    def __eq__(self, other):
        if type(self) != type(other): return False
        return self.ale_start == other.ale_start and self.ale_end == other.ale_end

    def get_start_stop(self):
        return[self.ale_start, self.ale_end]

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.symbol, self.gene, self.id, self.strand, "ALE", self.count,
                                                     self.ale_start, self.ale_end)
