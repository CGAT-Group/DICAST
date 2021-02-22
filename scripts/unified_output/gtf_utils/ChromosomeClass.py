import sys


class Chromosome:
    def __init__(self, id):
        self.id = id
        self.genes = {}

    def addGene(self, gene):
        if gene.getID() in self.genes:
            print("gene:" + gene.getID() + " already present in chromosome " + gene.getChromosome(), file=sys.stderr)
            sys.exit(1)
        self.genes[gene.getID()] = gene
