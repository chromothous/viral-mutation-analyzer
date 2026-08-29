import math

from classes.logger import Logger
from constants.motifs import MOTIFS


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
        canonical_bases = {"A", "T", "U", "G", "C"}
        known_counts = {
            nucleotide: count
            for nucleotide, count in counts.items()
            if nucleotide in canonical_bases
        }
        known_length = sum(known_counts.values())
        if known_length == 0:
            if "N" in counts:
                self.logger.warning("Sequence contains only ambiguous bases.")
                return 0.0
            raise ValueError("Sequence contains no recognized nucleotides.")
        entropy = 0.0
        for count in known_counts.values():
            probability = count / known_length
            entropy -= probability * math.log2(probability)
        maximum_entropy = math.log2(4)
        complexity = entropy / maximum_entropy
        complexity = min(max(complexity, 0.0), 1.0)
        if "N" in counts:
            self.logger.warning(
                f"Sequence complexity calculated with {counts['N']} ambiguous N bases excluded from known-base entropy."
            )
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

    def get_repeats(self, repeat_length):
        if self.sequence is None:
            raise ValueError("No sequence has been loaded.")
        if not isinstance(repeat_length, int):
            raise TypeError("Repeat length must be an integer.")
        if repeat_length <= 0:
            raise ValueError("Repeat length must be greater than zero.")
        if repeat_length > len(self.sequence):
            raise ValueError("Repeat length cannot exceed sequence length.")
        repeats = {}
        for start in range(len(self.sequence) - repeat_length + 1):
            pattern = self.sequence[start:start + repeat_length]
            repeats.setdefault(pattern, []).append(start + 1)
        repeated_patterns = {}
        for pattern, positions in repeats.items():
            if len(positions) > 1:
                repeated_patterns[pattern] = positions
        self.logger.info(
            f"Detected {len(repeated_patterns)} repeated patterns of length {repeat_length}."
        )
        return repeated_patterns

    def get_repeat_frequency(self, repeat_length):
        repeats = self.get_repeats(repeat_length)
        frequency = sum(len(positions) for positions in repeats.values())
        self.logger.info(
            f"Repeat frequency calculated: {frequency} occurrences."
        )
        return frequency

    def get_regional_repeat_density(self, window_size, repeat_length):
        regions = self.get_regions(window_size)
        regional_density = []
        for region in regions:
            sequence = region["sequence"]
            if len(sequence) < repeat_length:
                density = 0.0
            else:
                repeats = {}
                for start in range(len(sequence) - repeat_length + 1):
                    pattern = sequence[start:start + repeat_length]
                    repeats.setdefault(pattern, 0)
                    repeats[pattern] += 1
                repeated_occurrences = sum(
                    count for count in repeats.values() if count > 1
                )
                density = repeated_occurrences / len(sequence)
            regional_density.append({
                "region": region["region"],
                "start": region["start"],
                "end": region["end"],
                "repeat_density": density
            })
            if density > 0.5:
                self.logger.warning(
                    f"High repeat density detected: {region['start']}-{region['end']}."
                )
        self.logger.info(
            f"Regional repeat density calculated for {len(regional_density)} regions."
        )
        return regional_density

    def get_kmer_counts(self, k, sequence=None):
        if sequence is None:
            if self.sequence is None:
                raise ValueError("No sequence has been loaded.")
            sequence = self.sequence
        if not isinstance(k, int):
            raise TypeError("K-mer length must be an integer.")
        if k <= 0:
            raise ValueError("K-mer length must be greater than zero.")
        if k > len(sequence):
            raise ValueError("K-mer length cannot exceed sequence length.")
        counts = {}
        for start in range(len(sequence) - k + 1):
            kmer = sequence[start:start + k]
            if "N" in kmer:
                continue
            counts[kmer] = counts.get(kmer, 0) + 1
        self.logger.info(
            f"Calculated {len(counts)} unique k-mers using k={k}."
        )
        if "N" in sequence:
            self.logger.warning("K-mers containing ambiguous N bases were excluded.")
        return counts

    def get_kmer_frequencies(self, k, sequence=None):
        counts = self.get_kmer_counts(k, sequence)
        if sequence is None:
            sequence = self.sequence
        valid_kmers = sum(counts.values())
        if valid_kmers == 0:
            self.logger.warning("No valid k-mers were available for frequency analysis.")
            return {}
        frequencies = {
            kmer: count / valid_kmers
            for kmer, count in counts.items()
        }
        self.logger.info(
            f"K-mer frequencies calculated using k={k}."
        )
        return frequencies

    def get_kmer_diversity(self, k, sequence=None):
        counts = self.get_kmer_counts(k, sequence)
        if sequence is None:
            sequence = self.sequence
        valid_kmers = sum(counts.values())
        if valid_kmers == 0:
            self.logger.warning("No valid k-mers were available for diversity analysis.")
            return 0.0
        diversity = len(counts) / valid_kmers
        self.logger.info(
            f"K-mer diversity calculated using k={k}: {diversity:.4f}."
        )
        return diversity

    def get_regional_kmer_profiles(self, window_size, k):
        regions = self.get_regions(window_size)
        regional_profiles = []
        for region in regions:
            counts = self.get_kmer_counts(k, region["sequence"])
            total = sum(counts.values())
            frequencies = {}
            if total > 0:
                frequencies = {
                    kmer: count / total
                    for kmer, count in counts.items()
                }
            regional_profiles.append({
                "region": region["region"],
                "start": region["start"],
                "end": region["end"],
                "k": k,
                "counts": counts,
                "frequencies": frequencies,
                "unique_kmers": len(counts),
                "total_kmers": total
            })
        self.logger.info(
            f"Regional k-mer profiles calculated for {len(regional_profiles)} regions."
        )
        return regional_profiles

    def get_motif_positions(self, motif):
        if self.sequence is None:
            raise ValueError("No sequence has been loaded.")
        if not isinstance(motif, str):
            raise TypeError("Motif must be a string.")
        motif = motif.upper()
        if not motif:
            raise ValueError("Motif cannot be empty.")
        positions = []
        for start in range(len(self.sequence) - len(motif) + 1):
            if self.sequence[start:start + len(motif)] == motif:
                positions.append(start + 1)
        self.logger.info(
            f"Detected {len(positions)} occurrences of motif {motif}."
        )
        return positions

    def get_motif_analysis(self, motifs=None):
        if self.sequence is None:
            raise ValueError("No sequence has been loaded.")
        if motifs is None:
            motifs = MOTIFS
        if not isinstance(motifs, dict):
            raise TypeError("Motifs must be provided as a dictionary.")
        analysis = {}
        for name, motif in motifs.items():
            positions = self.get_motif_positions(motif)
            analysis[name] = {
                "motif": motif,
                "count": len(positions),
                "positions": positions,
                "density": len(positions) / self.get_length()
            }
        self.logger.info(
            f"Motif analysis completed for {len(analysis)} motifs."
        )
        return analysis

    def get_regional_motif_profiles(self, window_size, motifs=None):
        if motifs is None:
            motifs = MOTIFS
        if not isinstance(motifs, dict):
            raise TypeError("Motifs must be provided as a dictionary.")
        regions = self.get_regions(window_size)
        regional_profiles = []
        for region in regions:
            profile = {
                "region": region["region"],
                "start": region["start"],
                "end": region["end"],
                "motifs": {}
            }
            for name, motif in motifs.items():
                motif = motif.upper()
                positions = []
                for start in range(
                    region["start"] - 1,
                    min(region["end"], len(self.sequence) - len(motif) + 1)
                ):
                    if self.sequence[start:start + len(motif)] == motif:
                        positions.append(start + 1)
                profile["motifs"][name] = {
                    "motif": motif,
                    "count": len(positions),
                    "positions": positions,
                    "density": len(positions) / len(self.sequence)
                }
            regional_profiles.append(profile)
        self.logger.info(
            f"Regional motif profiles calculated for {len(regional_profiles)} regions."
        )
        return regional_profiles