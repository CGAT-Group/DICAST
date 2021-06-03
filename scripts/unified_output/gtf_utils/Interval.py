class Interval:

    def __init__(self, start, end):
        self.start = int(start)
        self.end = int(end)

    def __str__(self):
        return "[" + str(self.start) + ", " + str(self.end) + "]"

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __ne__(self, other):
        return not self.__eq__(other)

    def coveredBy(self, other):
        return self.start > other and self.end < other.end

    # also allows for swapped start and end
    def pseudo_equal(self, other):
        return self.__eq__(other) or (self.end == other.start and self.start == other.end)
