class Event:

    id = None
    start = 0
    end = 0
    gene = None
    symbol = None
    exon_list = None

    def __init__(self, idx, gene, symbol):
        self.id = idx
        self.gene = gene
        self.symbol = symbol

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
        super().__init__(idx, gene, symbol)
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
        super().__init__(idx, gene, symbol)
        self.mee_exons = mee_exons
        self.strand = strand
        self.symbol = symbol
        self.count = count

    def to_string(self):
        print("MEE-Event: " + str(self.id) + " mut. exclusive exons are: " + str(self.mee_exons))

    def __str__(self):
        mee_exons_starts = ""
        mee_exons_ends = ""
        for i in range(len(self.mee_exons)-1):
            mee_exons_starts += str(self.mee_exons[i][0] + ",")
            mee_exons_ends += str(self.mee_exons[i][1] + ",")
        mee_exons_starts += str(self.mee_exons[len(self.mee_exons)-1][0])
        mee_exons_ends += str(self.mee_exons[len(self.mee_exons)-1][1])
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.symbol, self.gene, self.id, self.strand, "MEE", self.count, mee_exons_starts, mee_exons_ends)

    def __hash__(self):
        return hash((self.gene, frozenset(frozenset(ev) for ev in self.mee_exons)))

    def __eq__(self, other):
        sorted_exons = sorted([list(map(int, i)) for i in self.mee_exons])
        other_sorted_exons = sorted([list(map(int, i)) for i in other.mee_exons])
        if type(self) != type(other): return False
        if len(sorted_exons) != len(other_sorted_exons): return False
        b = True
        for i in range(len(sorted_exons)):
            if sorted_exons[i] != other_sorted_exons[i]:
                b = False
        return b

    def get_start_stop(self):
        return self.mee_exons


class MesEvent(Event):

    mes_skipped_exons = []
    strand = ""
    count = 0

    def __init__(self, idx,  mes_skipped, strand, gene, symbol, count=1):
        super().__init__(idx, gene, symbol)
        self.mes_skipped_exons = mes_skipped
        self.strand = strand
        self.symbol = symbol
        self.count = count

    def __eq__(self, other):
        sorted_exons = sorted([list(map(int, i)) for i in self.mes_skipped_exons])
        other_sorted_exons = sorted([list(map(int, i)) for i in other.mes_skipped_exons])
        if type(self) != type(other): return False
        if len(sorted_exons) != len(other_sorted_exons): return False
        b = True
        for i in range(len(sorted_exons)):
            if sorted_exons[i] != other_sorted_exons[i]:
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
        super().__init__(idx, gene, symbol)
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
        super().__init__(idx, gene, symbol)
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
        super().__init__(idx, gene, symbol)
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
        super().__init__(idx, gene, symbol)
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
        super().__init__(idx, gene, symbol)
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