from classes.logger import Logger


class Analyzer:
    DNA_BASES = set("ATGCN")
    RNA_BASES = set("AUGCN")

    def __init__(self, sequence=None):
        self.logger = Logger()
        self.sequence = None
        self.sequence_type = None
        self.logger.info("Analyzer initialized.")
        if sequence is not None:
            self.set_sequence(sequence)

    def set_sequence(self, sequence):
        if not isinstance(sequence, str):
            self.logger.failure("Sequence rejected: input is not a string.")
            raise TypeError("Sequence must be a string.")
        sequence = sequence.upper().replace(" ", "").replace("\n", "").replace("\r", "")
        if not sequence:
            self.logger.failure("Sequence rejected: sequence is empty.")
            raise ValueError("Sequence cannot be empty.")
        dna_valid = set(sequence).issubset(self.DNA_BASES)
        rna_valid = set(sequence).issubset(self.RNA_BASES)
        if dna_valid and "T" in sequence:
            self.sequence_type = "DNA"
        elif rna_valid and "U" in sequence:
            self.sequence_type = "RNA"
        elif dna_valid:
            self.sequence_type = "DNA"
            self.logger.warning("Sequence contains only ambiguous/resolved DNA-compatible bases.")
        else:
            self.logger.failure("Sequence rejected: invalid nucleotide detected.")
            raise ValueError("Sequence contains invalid nucleotides.")
        self.sequence = sequence
        self.logger.info(f"{self.sequence_type} sequence accepted ({len(sequence)} bases).")

    def get_sequence(self):
        return self.sequence

    def get_sequence_type(self):
        return self.sequence_type