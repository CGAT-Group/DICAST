from tool_parser.SGSeqReader import SGSeqReader
from gtf_utils.GTFParser import GTFParser


class SGSeqAnnoReader(SGSeqReader):

    NAME = "SGSeq_Anno"

    def __init__(self, filepath, gtf: GTFParser, combine_me):
        super().__init__(filepath, gtf, combine_me)

