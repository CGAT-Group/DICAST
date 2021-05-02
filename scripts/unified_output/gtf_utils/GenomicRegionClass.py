class GenomicRegion:

    def __init__(self, feature_id, start, end, on_positive_strand):
        self.feature_id = feature_id
        self.start = start
        self.end = end
        self.on_positive_strand = on_positive_strand

    def getID(self):
        return self.feature_id

    def getStart(self):
        return self.start

    def getEnd(self):
        return self.end

    def is_on_positive_strand(self):
        return self.on_positive_strand

    def getLength(self):
        return self.end - self.start + 1

    def __str__(self):
        return self.start + "-" + self.end
