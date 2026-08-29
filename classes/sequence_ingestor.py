from classes.logger import Logger


class SequenceIngestor:
    DNA_BASES = set("ATGCN")
    RNA_BASES = set("AUGCN")

    def __init__(self):
        self.logger = Logger()
        self.sequence = None
        self.sequence_type = None

    def ingest(self, sequence):
        if not isinstance(sequence, str):
            self.logger.failure("Sequence ingestion failed: input is not a string.")
            raise TypeError("Sequence must be a string.")

        sequence = sequence.upper()
        sequence = "".join(sequence.split())

        if not sequence:
            self.logger.failure("Sequence ingestion failed: sequence is empty.")
            raise ValueError("Sequence cannot be empty.")

        dna_valid = set(sequence).issubset(self.DNA_BASES)
        rna_valid = set(sequence).issubset(self.RNA_BASES)

        if dna_valid and "T" in sequence:
            self.sequence_type = "DNA"
        elif rna_valid and "U" in sequence:
            self.sequence_type = "RNA"
        elif dna_valid and rna_valid:
            self.sequence_type = "UNKNOWN"
            self.logger.warning(
                "Sequence contains only bases shared by DNA and RNA."
            )
        else:
            self.logger.failure(
                "Sequence ingestion failed: invalid nucleotide detected."
            )
            raise ValueError("Sequence contains invalid nucleotides.")

        self.sequence = sequence

        self.logger.info(
            f"Sequence ingestion successful: {self.sequence_type} "
            f"sequence containing {len(self.sequence)} bases."
        )

        if "N" in self.sequence:
            self.logger.warning(
                f"Sequence contains {self.sequence.count('N')} ambiguous N bases."
            )

        return {
            "sequence": self.sequence,
            "sequence_type": self.sequence_type,
            "length": len(self.sequence),
            "ambiguous_bases": self.sequence.count("N")
        }

    def get_sequence(self):
        return self.sequence

    def get_sequence_type(self):
        return self.sequence_type