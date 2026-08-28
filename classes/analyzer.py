class Analyzer:
    DNA_BASES = set("ATGCN")
    RNA_BASES = set("AUGCN")

    def __init__(self, sequence=None):
        self.sequence = None
        self.sequence_type = None

        if sequence is not None:
            self.set_sequence(sequence)

    def set_sequence(self, sequence):
        if not isinstance(sequence, str):
            raise TypeError("Sequence must be a string.")

        sequence = sequence.upper().replace(" ", "").replace("\n", "").replace("\r", "")

        if not sequence:
            raise ValueError("Sequence cannot be empty.")

        dna_valid = set(sequence).issubset(self.DNA_BASES)
        rna_valid = set(sequence).issubset(self.RNA_BASES)

        if dna_valid and "T" in sequence:
            self.sequence_type = "DNA"
        elif rna_valid and "U" in sequence:
            self.sequence_type = "RNA"
        elif dna_valid:
            self.sequence_type = "DNA"
        else:
            raise ValueError("Sequence contains invalid nucleotides.")

        self.sequence = sequence

    def get_sequence(self):
        return self.sequence

    def get_sequence_type(self):
        return self.sequence_type