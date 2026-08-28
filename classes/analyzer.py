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

    def get_length(self):
        if self.sequence is None:
            raise ValueError("No sequence has been loaded.")
        return len(self.sequence)

    def get_nucleotide_counts(self):
        if self.sequence is None:
            raise ValueError("No sequence has been loaded.")
        counts = {
            "A": self.sequence.count("A"),
            "T": self.sequence.count("T"),
            "U": self.sequence.count("U"),
            "G": self.sequence.count("G"),
            "C": self.sequence.count("C"),
            "N": self.sequence.count("N")
        }
        self.logger.info("Nucleotide counts calculated.")
        if counts["N"] > 0:
            self.logger.warning(f"Sequence contains {counts['N']} ambiguous N bases.")
        return counts

    def get_nucleotide_frequencies(self):
        counts = self.get_nucleotide_counts()
        length = self.get_length()
        frequencies = {}
        for nucleotide, count in counts.items():
            frequencies[nucleotide] = count / length
        self.logger.info("Nucleotide frequencies calculated.")
        return frequencies

    def get_gc_content(self):
        counts = self.get_nucleotide_counts()
        length = self.get_length()
        gc_content = (counts["G"] + counts["C"]) / length
        self.logger.info(f"GC content calculated: {gc_content:.4f}.")
        return gc_content