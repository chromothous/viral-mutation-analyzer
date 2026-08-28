import math

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

    def get_regions(self, window_size):
        if self.sequence is None:
            raise ValueError("No sequence has been loaded.")
        if not isinstance(window_size, int):
            raise TypeError("Window size must be an integer.")
        if window_size <= 0:
            raise ValueError("Window size must be greater than zero.")
        regions = []
        for start in range(0, len(self.sequence), window_size):
            end = min(start + window_size, len(self.sequence))
            regions.append({
                "region": len(regions) + 1,
                "start": start + 1,
                "end": end,
                "sequence": self.sequence[start:end]
            })
        self.logger.info(
            f"Sequence divided into {len(regions)} regions using a window size of {window_size}."
        )
        return regions

    def get_regional_statistics(self, window_size):
        regions = self.get_regions(window_size)
        regional_statistics = []
        for region in regions:
            sequence = region["sequence"]
            length = len(sequence)
            counts = {
                "A": sequence.count("A"),
                "T": sequence.count("T"),
                "U": sequence.count("U"),
                "G": sequence.count("G"),
                "C": sequence.count("C"),
                "N": sequence.count("N")
            }
            frequencies = {
                nucleotide: count / length
                for nucleotide, count in counts.items()
            }
            gc_content = (counts["G"] + counts["C"]) / length
            regional_statistics.append({
                "region": region["region"],
                "start": region["start"],
                "end": region["end"],
                "length": length,
                "counts": counts,
                "frequencies": frequencies,
                "gc_content": gc_content,
                "ambiguous_bases": counts["N"]
            })
        self.logger.info(
            f"Regional statistics calculated for {len(regional_statistics)} regions."
        )
        return regional_statistics

    def get_sequence_complexity(self, sequence=None):
        if sequence is None:
            if self.sequence is None:
                raise ValueError("No sequence has been loaded.")
            sequence = self.sequence
        if not isinstance(sequence, str):
            raise TypeError("Sequence must be a string.")
        if not sequence:
            raise ValueError("Sequence cannot be empty.")
        counts = {}
        for nucleotide in sequence:
            counts[nucleotide] = counts.get(nucleotide, 0) + 1
        length = len(sequence)
        entropy = 0.0
        for count in counts.values():
            probability = count / length
            entropy -= probability * math.log2(probability)
        maximum_entropy = math.log2(4)
        complexity = entropy / maximum_entropy
        self.logger.info(f"Sequence complexity calculated: {complexity:.4f}.")
        return complexity

    def get_regional_complexity(self, window_size):
        regions = self.get_regions(window_size)
        regional_complexity = []
        for region in regions:
            complexity = self.get_sequence_complexity(region["sequence"])
            regional_complexity.append({
                "region": region["region"],
                "start": region["start"],
                "end": region["end"],
                "complexity": complexity
            })
            if complexity < 0.5:
                self.logger.warning(
                    f"Low-complexity region detected: {region['start']}-{region['end']}."
                )
        self.logger.info(
            f"Regional complexity calculated for {len(regional_complexity)} regions."
        )
        return regional_complexity