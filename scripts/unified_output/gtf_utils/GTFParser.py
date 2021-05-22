import gzip

from gtf_utils.ChromosomeClass import Chromosome
from gtf_utils.ExonClass import Exon
from gtf_utils.GeneClass import Gene
from gtf_utils.TranscriptClass import Transcript


class GTFParser:
    GENE_IDENTIFIER = "gene"
    TRANSCRIPT_IDENTIFIER = "transcript"
    EXON_IDENTIFIER = "exon"

    def __init__(self, gtf_filepath):

        self.chromosomes = {}
        self.genes = {}
        self.transcripts = {}

        current_gene = None
        current_transcript = None

        open_fct = gzip.open if gtf_filepath.endswith('.gz') else open

        with open_fct(gtf_filepath) as gtf:
            for line in gtf:
                if line.startswith("#"):
                    continue

                # remove newline characters from line and split by tab
                # remove trailing ";" used in last column of gene lines in gtf
                line = line.strip().rstrip(";").split(sep="\t")

                chr_id = line[0]
                if chr_id in self.chromosomes:
                    current_chromosome = self.chromosomes[chr_id]
                else:
                    current_chromosome = Chromosome(chr_id)
                    self.chromosomes[chr_id] = current_chromosome

                feature_type = line[2]

                if feature_type == GTFParser.EXON_IDENTIFIER:
                    # gene, transcript, exon_number, start, end, tr_start, tr_end, is_on_positive_strand,
                    # is_template=True
                    exon_attributes = line[8]
                    exon_attributes = parseAttributes(exon_attributes)
                    is_template = exon_attributes["template"] == "TRUE"
                    current_exon = Exon(current_gene, current_transcript,
                                        int(exon_attributes["gene_exon_number"]), int(line[3]), int(line[4]),
                                        int(exon_attributes["tr_start"]), int(exon_attributes["tr_end"]),
                                        line[6] == "+", is_template)
                    current_transcript.addExon(current_exon)
                    if is_template or not current_exon.getID() in current_gene.exons:
                        current_gene.addExon(current_exon)

                elif feature_type == GTFParser.TRANSCRIPT_IDENTIFIER:
                    # id, chr, gene, start, end, is_on_positive_strand, tr_start, tr_end, is_template = True
                    transcript_attributes = line[8]
                    transcript_attributes = parseAttributes(transcript_attributes)

                    current_transcript = Transcript(transcript_attributes["transcript_id"], current_chromosome,
                                                    current_gene, int(line[3]), int(line[4]),
                                                    line[6] == "+",
                                                    int(transcript_attributes["tr_start"]),
                                                    int(transcript_attributes["tr_end"]),
                                                    transcript_attributes["template"] == "TRUE")

                    if not current_transcript.is_template:
                        event_type = transcript_attributes["transcript_id"].split("_")[1]
                        current_transcript.addEventType(event_type)

                    current_gene.addTranscript(current_transcript)
                    self.transcripts[current_transcript.getID()] = current_transcript

                elif feature_type == GTFParser.GENE_IDENTIFIER:
                    # id, chr, start, end, is_on_positive_strand
                    gene_id = line[8].split()[1].replace("\"", "")
                    current_gene = Gene(gene_id, current_chromosome, int(line[3]), int(line[4]), line[6] == "+")
                    current_chromosome.addGene(current_gene)
                    self.genes[gene_id] = current_gene

        # sort transcripts by position not number
        for tr in self.transcripts.values():
            if not tr.is_on_positive_strand():
                tr.getExons().reverse()

        for g in self.genes.values():
            g.calcEventJunctions()


# return dict from features attributes (input gtf_line[8])
def parseAttributes(feature_attributes):
    attributes = {}
    feature_attributes = feature_attributes.replace("\"", "").split("; ")
    for attr in feature_attributes:
        attr = attr.split()
        attributes[attr[0]] = attr[1]

    return attributes


def parseEventType(transcript_id):
    return GTFParser.EVENT_TYPES[transcript_id.split("_")[1]]
