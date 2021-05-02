from gtf_utils.GenomicRegionClass import GenomicRegion

# exons sorted by genomic position when read via GTFParser!!!
class Transcript(GenomicRegion):
    def __init__(self, feature_id, chr, gene, start, end, is_on_positive_strand, tr_start, tr_end, is_template=True):
        super().__init__(feature_id, start, end, is_on_positive_strand)
        self.chr = chr
        self.gene = gene
        self.tr_start = tr_start
        self.tr_end = tr_end
        self.is_template = is_template
        self.exons = list()
        self.event_types = list()

    def addExon(self, exon):
        self.exons.append(exon)

    def getChromosome(self):
        return self.chr

    def getExons(self):
        return self.exons

    def getExon(self, exon_number):
        for e in self.exons:
            if e.getID() == exon_number:
                return e
        return None

    def addEventType(self, event_type):
        self.event_types.append(event_type)

    def getEventTypes(self):
        return self.event_types
