from gtf_utils.GenomicRegionClass import GenomicRegion


class Exon(GenomicRegion):

    def __init__(self, gene, transcript, exon_number, start, end, tr_start, tr_end, is_on_positive_strand,
                 is_template=True, is_artificial=False):
        super().__init__(exon_number, start, end, is_on_positive_strand)
        self.tr_start = tr_start
        self.tr_end = tr_end
        self.is_on_positive_strand = is_on_positive_strand
        self.is_template = is_template
        self.is_artificial = is_artificial

    def __str__(self):
        return str(self.start) + "-" + str(self.end)

    def equals(self, start, end):
        if int(start) == self.start and int(end) == self.end:
            return True
        return False

    def get_artificial_coords(self):
        if self.start == "nan":
            return int(self.end)-1, int(self.end)
        else:
            return int(self.start), int(self.start)+1

    def set_start(self, start):
        self.start = start

    def set_end(self, end):
        self.end = end