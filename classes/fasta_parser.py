from classes.logger import Logger
from classes.sequence_ingestor import SequenceIngestor


class FastaParser:
    def __init__(self):
        self.logger = Logger()
        self.ingestor = SequenceIngestor()
        self.header = None
        self.sequence = None
        self.metadata = None

    def parse(self, content):
        if not isinstance(content, str):
            self.logger.failure("FASTA parsing failed: input is not a string.")
            raise TypeError("FASTA content must be a string.")

        lines = content.strip().splitlines()

        if not lines:
            self.logger.failure("FASTA parsing failed: input is empty.")
            raise ValueError("FASTA content cannot be empty.")

        if not lines[0].startswith(">"):
            self.logger.failure("FASTA parsing failed: missing header.")
            raise ValueError("FASTA content must begin with a header.")

        header = lines[0][1:].strip()

        if not header:
            self.logger.failure("FASTA parsing failed: header is empty.")
            raise ValueError("FASTA header cannot be empty.")

        sequence_lines = lines[1:]

        if not sequence_lines:
            self.logger.failure("FASTA parsing failed: sequence is missing.")
            raise ValueError("FASTA sequence cannot be empty.")

        sequence = "".join(sequence_lines)

        result = self.ingestor.ingest(sequence)

        self.header = header
        self.sequence = result["sequence"]
        self.metadata = {
            "header": self.header,
            "sequence": self.sequence,
            "sequence_type": result["sequence_type"],
            "length": result["length"],
            "ambiguous_bases": result["ambiguous_bases"]
        }

        self.logger.info(
            f"FASTA record parsed successfully: {self.header}."
        )

        return self.metadata

    def get_header(self):
        return self.header

    def get_sequence(self):
        return self.sequence

    def get_metadata(self):
        return self.metadata