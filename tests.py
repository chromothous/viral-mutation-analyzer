import os
import sys
import traceback

from classes.analyzer import Analyzer
from classes.logger import Logger
from classes.sequence_ingestor import SequenceIngestor
from classes.fasta_parser import FastaParser
from classes.configuration import Configuration
from classes.feature_normalizer import FeatureNormalizer

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def green(text):
    return f"\033[32m{text}\033[0m"

def red(text):
    return f"\033[31m{text}\033[0m"

def full_test():
    tests = 0
    success = 0
    failure = 0

    try:
        tests += 1
        assert os.path.isdir("classes")
        assert os.path.isdir("constants")
        assert os.path.isfile("main.py")
        assert os.path.isfile("classes/analyzer.py")
        assert os.path.isfile("tests.py")
        analyzer = Analyzer()
        assert isinstance(analyzer, Analyzer)
        print(green("Version 0.0.1 foundation is online."))
        success += 1
    except Exception as e:
        failure += 1
        red(traceback.print_exc())
        print(red(e))
        print(red("Version 0.0.1 foundation failed."))

    try:
        tests += 1
        analyzer = Analyzer("ATGCNNATGC")
        assert analyzer.get_sequence() == "ATGCNNATGC"
        assert analyzer.get_sequence_type() == "DNA"
        analyzer = Analyzer("AUGCNN AUGC")
        assert analyzer.get_sequence() == "AUGCNNAUGC"
        assert analyzer.get_sequence_type() == "RNA"
        analyzer.set_sequence("ATGCNN")
        assert analyzer.get_sequence() == "ATGCNN"
        assert analyzer.get_sequence_type() == "DNA"
        try:
            Analyzer("ATGCX")
            raise AssertionError("Invalid nucleotide was accepted.")
        except ValueError:
            pass
        try:
            Analyzer("")
            raise AssertionError("Empty sequence was accepted.")
        except ValueError:
            pass
        try:
            Analyzer(12345)
            raise AssertionError("Non-string sequence was accepted.")
        except TypeError:
            pass
        print(green("Version 0.0.2 sequence foundation is online."))
        success += 1
    except Exception as e:
        failure += 1
        red(traceback.print_exc())
        print(red(e))
        print(red("Version 0.0.2 sequence foundation failed."))

    try:
        tests += 1
        logger = Logger()
        assert hasattr(logger, "info")
        assert hasattr(logger, "warning")
        assert hasattr(logger, "failure")
        analyzer = Analyzer("ATGCNNATGC")
        assert isinstance(analyzer.logger, Logger)
        assert analyzer.get_sequence() == "ATGCNNATGC"
        assert analyzer.get_sequence_type() == "DNA"
        analyzer.logger.info("Testing logger information output.")
        analyzer.logger.warning("Testing logger warning output.")
        analyzer.logger.failure("Testing logger failure output.")
        print(green("Version 0.0.3 logging foundation is online."))
        success += 1
    except Exception as e:
        failure += 1
        red(traceback.print_exc())
        print(red(e))
        print(red("Version 0.0.3 logging foundation failed."))

    try:
        tests += 1
        analyzer = Analyzer("ATGCCNGTAA")
        assert analyzer.get_length() == 10
        counts = analyzer.get_nucleotide_counts()
        assert counts["A"] == 3
        assert counts["T"] == 2
        assert counts["G"] == 2
        assert counts["C"] == 2
        assert counts["N"] == 1
        assert counts["U"] == 0
        frequencies = analyzer.get_nucleotide_frequencies()
        assert frequencies["A"] == 0.3
        assert frequencies["T"] == 0.2
        assert frequencies["G"] == 0.2
        assert frequencies["C"] == 0.2
        assert frequencies["N"] == 0.1
        assert frequencies["U"] == 0.0
        assert analyzer.get_gc_content() == 0.4
        print(green("Version 0.0.4 sequence statistics are online."))
        success += 1
    except Exception as e:
        failure += 1
        red(traceback.print_exc())
        print(red(e))
        print(red("Version 0.0.4 sequence statistics failed."))

    try:
        tests += 1
        analyzer = Analyzer("ATGCNNATGCATGC")
        regions = analyzer.get_regions(5)
        assert isinstance(regions, list)
        assert len(regions) == 3
        assert regions[0]["region"] == 1
        assert regions[0]["start"] == 1
        assert regions[0]["end"] == 5
        assert regions[0]["sequence"] == "ATGCN"
        assert regions[1]["region"] == 2
        assert regions[1]["start"] == 6
        assert regions[1]["end"] == 10
        assert regions[1]["sequence"] == "NATGC"
        assert regions[2]["region"] == 3
        assert regions[2]["start"] == 11
        assert regions[2]["end"] == 14
        assert regions[2]["sequence"] == "ATGC"
        statistics = analyzer.get_regional_statistics(5)
        assert isinstance(statistics, list)
        assert len(statistics) == 3
        assert statistics[0]["length"] == 5
        assert statistics[0]["counts"]["N"] == 1
        assert statistics[1]["length"] == 5
        assert statistics[1]["counts"]["N"] == 1
        assert statistics[2]["length"] == 4
        assert statistics[2]["counts"]["A"] == 1
        assert statistics[2]["counts"]["T"] == 1
        assert statistics[2]["counts"]["G"] == 1
        assert statistics[2]["counts"]["C"] == 1
        assert statistics[2]["counts"]["N"] == 0
        assert statistics[2]["gc_content"] == 0.5
        try:
            analyzer.get_regions(0)
            raise AssertionError("Zero window size was accepted.")
        except ValueError:
            pass
        try:
            analyzer.get_regions(-5)
            raise AssertionError("Negative window size was accepted.")
        except ValueError:
            pass
        try:
            analyzer.get_regions("5")
            raise AssertionError("Non-integer window size was accepted.")
        except TypeError:
            pass
        print(green("Version 0.0.5 regional sequence analysis is online."))
        success += 1
    except Exception as e:
        failure += 1
        red(traceback.print_exc())
        print(red(e))
        print(red("Version 0.0.5 regional sequence analysis failed."))

    try:
        tests += 1
        analyzer = Analyzer("AAAACCCCGGGGTTTT")
        complexity = analyzer.get_sequence_complexity()
        assert isinstance(complexity, float), "Complexity is not a float."
        assert complexity == 1.0, "Four equally distributed bases should have complexity 1.0."
        analyzer = Analyzer("AAAAAAAAAAAAAAAA")
        complexity = analyzer.get_sequence_complexity()
        assert complexity == 0.0, "A single repeated base should have complexity 0.0."
        analyzer = Analyzer("ATGCATGCATGC")
        complexity = analyzer.get_sequence_complexity()
        assert complexity == 1.0, "Four equally distributed bases should have complexity 1.0."
        analyzer = Analyzer("AAAAATTTTT")
        complexity = analyzer.get_sequence_complexity()
        assert complexity == 0.5, "Two equally distributed bases should have complexity 0.5."
        regional_complexity = analyzer.get_regional_complexity(5)
        assert isinstance(regional_complexity, list), "Regional complexity did not return a list."
        assert len(regional_complexity) == 2, "Expected two complexity regions."
        assert regional_complexity[0]["region"] == 1, "First complexity region number is incorrect."
        assert regional_complexity[0]["start"] == 1, "First complexity region start is incorrect."
        assert regional_complexity[0]["end"] == 5, "First complexity region end is incorrect."
        assert regional_complexity[0]["complexity"] == 0.0, "First region complexity is incorrect."
        assert regional_complexity[1]["region"] == 2, "Second complexity region number is incorrect."
        assert regional_complexity[1]["start"] == 6, "Second complexity region start is incorrect."
        assert regional_complexity[1]["end"] == 10, "Second complexity region end is incorrect."
        assert regional_complexity[1]["complexity"] == 0.0, "Second region complexity is incorrect."
        analyzer = Analyzer("ATGCNN")
        complexity = analyzer.get_sequence_complexity()
        assert isinstance(complexity, float), "Complexity with N is not a float."
        assert 0.0 <= complexity <= 1.0, "Complexity with N is outside the expected range."
        print(green("Version 0.0.6 sequence complexity analysis is online."))
        success += 1
    except Exception as e:
        failure += 1
        print(red(e))
        print(red("Version 0.0.6 sequence complexity analysis failed."))

    try:
        tests += 1
        analyzer = Analyzer("ATGCATGCATGC")
        repeats = analyzer.get_repeats(4)
        assert isinstance(repeats, dict), "Repeat analysis did not return a dictionary."
        assert "ATGC" in repeats, "Expected ATGC repeat was not detected."
        assert len(repeats["ATGC"]) == 3, "Incorrect number of ATGC repeat occurrences."
        assert repeats["ATGC"] == [1, 5, 9], "ATGC repeat positions are incorrect."
        frequency = analyzer.get_repeat_frequency(4)
        assert frequency == 9, "Repeat frequency should count all repeated pattern occurrences."
        analyzer = Analyzer("AAAAATTTTT")
        repeats = analyzer.get_repeats(2)
        assert "AA" in repeats, "AA repeat was not detected."
        assert "TT" in repeats, "TT repeat was not detected."
        assert len(repeats["AA"]) == 4, "Incorrect AA repeat count."
        assert len(repeats["TT"]) == 4, "Incorrect TT repeat count."
        regional_density = analyzer.get_regional_repeat_density(5, 2)
        assert isinstance(regional_density, list), "Regional repeat density did not return a list."
        assert len(regional_density) == 2, "Incorrect number of regional repeat-density results."
        assert regional_density[0]["region"] == 1, "First repeat-density region is incorrect."
        assert regional_density[1]["region"] == 2, "Second repeat-density region is incorrect."
        assert regional_density[0]["repeat_density"] > 0.0, "First region repeat density should be greater than zero."
        assert regional_density[1]["repeat_density"] > 0.0, "Second region repeat density should be greater than zero."
        try:
            analyzer.get_repeats(0)
            raise AssertionError("Zero repeat length was accepted.")
        except ValueError:
            pass
        try:
            analyzer.get_repeats(-1)
            raise AssertionError("Negative repeat length was accepted.")
        except ValueError:
            pass
        try:
            analyzer.get_repeats("2")
            raise AssertionError("Non-integer repeat length was accepted.")
        except TypeError:
            pass
        print(green("Version 0.0.7 repetitive sequence analysis is online."))
        success += 1
    except Exception as e:
        failure += 1
        print(red(e))
        print(red("Version 0.0.7 repetitive sequence analysis failed."))

    try:
        tests += 1
        analyzer = Analyzer("ATGCATGC")
        counts = analyzer.get_kmer_counts(2)
        assert isinstance(counts, dict), "K-mer counts did not return a dictionary."
        assert counts["AT"] == 2, "AT k-mer count is incorrect."
        assert counts["TG"] == 2, "TG k-mer count is incorrect."
        assert counts["GC"] == 2, "GC k-mer count is incorrect."
        assert counts["CA"] == 1, "CA k-mer count is incorrect."
        assert len(counts) == 4, "Unique k-mer count is incorrect."
        frequencies = analyzer.get_kmer_frequencies(2)
        assert abs(frequencies["AT"] - (2 / 7)) < 0.000001, "AT k-mer frequency is incorrect."
        assert abs(frequencies["TG"] - (2 / 7)) < 0.000001, "TG k-mer frequency is incorrect."
        assert abs(frequencies["GC"] - (2 / 7)) < 0.000001, "GC k-mer frequency is incorrect."
        assert abs(frequencies["CA"] - (1 / 7)) < 0.000001, "CA k-mer frequency is incorrect."
        diversity = analyzer.get_kmer_diversity(2)
        assert abs(diversity - (4 / 7)) < 0.000001, "K-mer diversity is incorrect."
        analyzer = Analyzer("ATGCNNATGC")
        counts = analyzer.get_kmer_counts(2)
        assert "NN" not in counts, "K-mer containing N was included."
        assert "NA" not in counts, "K-mer containing N was included."
        assert "CN" not in counts, "K-mer containing N was included."
        assert counts["AT"] == 2, "AT count with ambiguous bases is incorrect."
        regional_profiles = analyzer.get_regional_kmer_profiles(5, 2)
        assert isinstance(regional_profiles, list), "Regional k-mer profiles did not return a list."
        assert len(regional_profiles) == 2, "Incorrect number of regional k-mer profiles."
        assert regional_profiles[0]["region"] == 1, "First k-mer region number is incorrect."
        assert regional_profiles[1]["region"] == 2, "Second k-mer region number is incorrect."
        assert regional_profiles[0]["k"] == 2, "First region k value is incorrect."
        assert regional_profiles[1]["k"] == 2, "Second region k value is incorrect."
        try:
            analyzer.get_kmer_counts(0)
            raise AssertionError("Zero k-mer length was accepted.")
        except ValueError:
            pass
        try:
            analyzer.get_kmer_counts(-1)
            raise AssertionError("Negative k-mer length was accepted.")
        except ValueError:
            pass
        try:
            analyzer.get_kmer_counts("2")
            raise AssertionError("Non-integer k-mer length was accepted.")
        except TypeError:
            pass
        print(green("Version 0.0.8 k-mer analysis is online."))
        success += 1
    except Exception as e:
        failure += 1
        print(red(e))
        print(red("Version 0.0.8 k-mer analysis failed."))

    try:
        tests += 1
        analyzer = Analyzer("ATGCGCATGCGCATGC")
        positions = analyzer.get_motif_positions("ATGC")
        assert isinstance(positions, list), "Motif positions did not return a list."
        assert positions == [1, 7, 13], "Motif positions are incorrect."
        analysis = analyzer.get_motif_analysis({
            "START": "ATG",
            "GC": "GC"
        })
        assert isinstance(analysis, dict), "Motif analysis did not return a dictionary."
        assert "START" in analysis, "START motif is missing from analysis."
        assert "GC" in analysis, "GC motif is missing from analysis."
        assert analysis["START"]["motif"] == "ATG", "START motif value is incorrect."
        assert analysis["START"]["count"] == 3, "START motif count is incorrect."
        assert analysis["START"]["positions"] == [1, 7, 13], "START motif positions are incorrect."
        assert analysis["GC"]["count"] == 5, "GC motif count is incorrect."
        assert analysis["GC"]["positions"] == [3, 5, 9, 11, 15], "GC motif positions are incorrect."      
        assert abs(analysis["START"]["density"] - (3 / 16)) < 0.000001, "START motif density is incorrect."
        regional_profiles = analyzer.get_regional_motif_profiles(8, {
            "START": "ATG",
            "GC": "GC"
        })
        assert isinstance(regional_profiles, list), "Regional motif profiles did not return a list."
        assert len(regional_profiles) == 2, "Incorrect number of regional motif profiles."
        assert regional_profiles[0]["region"] == 1, "First motif region number is incorrect."
        assert regional_profiles[1]["region"] == 2, "Second motif region number is incorrect."
        assert regional_profiles[0]["motifs"]["START"]["count"] == 2, "First region START motif count is incorrect."
        assert regional_profiles[1]["motifs"]["START"]["count"] == 1, "Second region START motif count is incorrect."
        analyzer = Analyzer("ATGCNNATGC")
        positions = analyzer.get_motif_positions("ATGC")
        assert positions == [1, 7], "Motif detection with N is incorrect."
        try:
            analyzer.get_motif_positions("")
            raise AssertionError("Empty motif was accepted.")
        except ValueError:
            pass
        try:
            analyzer.get_motif_positions(123)
            raise AssertionError("Non-string motif was accepted.")
        except TypeError:
            pass
        print(green("Version 0.0.9 motif analysis is online."))
        success += 1
    except Exception as e:
        failure += 1
        print(red(e))
        print(red("Version 0.0.9 motif analysis failed."))

    try:
        tests += 1
        analyzer = Analyzer("ATGCATGCNNATGCATGC")
        analysis = analyzer.get_unified_regional_analysis(
            6,
            2,
            2,
            {
                "ATG": "ATG",
                "GC": "GC"
            }
        )
        assert isinstance(analysis, list), "Unified analysis did not return a list."
        assert len(analysis) == 3, "Unified analysis returned an incorrect number of regions."
        assert analysis[0]["region"] == 1, "First unified region number is incorrect."
        assert analysis[1]["region"] == 2, "Second unified region number is incorrect."
        assert analysis[2]["region"] == 3, "Third unified region number is incorrect."
        assert analysis[0]["start"] == 1, "First unified region start is incorrect."
        assert analysis[0]["end"] == 6, "First unified region end is incorrect."
        assert analysis[1]["start"] == 7, "Second unified region start is incorrect."
        assert analysis[1]["end"] == 12, "Second unified region end is incorrect."
        assert analysis[2]["start"] == 13, "Third unified region start is incorrect."
        assert analysis[2]["end"] == 18, "Third unified region end is incorrect."
        assert "statistics" in analysis[0], "Regional statistics are missing from unified analysis."
        assert "complexity" in analysis[0], "Regional complexity is missing from unified analysis."
        assert "repeat_density" in analysis[0], "Repeat density is missing from unified analysis."
        assert "kmer_profile" in analysis[0], "K-mer profile is missing from unified analysis."
        assert "motif_profile" in analysis[0], "Motif profile is missing from unified analysis."
        assert analysis[0]["statistics"]["length"] == 6, "First region statistics length is incorrect."
        assert analysis[1]["statistics"]["counts"]["N"] == 2, "Second region N count is incorrect."
        assert 0.0 <= analysis[0]["complexity"]["complexity"] <= 1.0, "First region complexity is outside the expected range."
        assert 0.0 <= analysis[1]["complexity"]["complexity"] <= 1.0, "Second region complexity is outside the expected range."
        assert analysis[0]["kmer_profile"]["k"] == 2, "First region k-mer length is incorrect."
        assert analysis[1]["kmer_profile"]["k"] == 2, "Second region k-mer length is incorrect."
        assert "ATG" in analysis[0]["motif_profile"]["motifs"], "ATG motif is missing from first region."
        assert "GC" in analysis[0]["motif_profile"]["motifs"], "GC motif is missing from first region."
        assert isinstance(analysis[0]["repeat_density"]["repeat_density"], float), "Repeat density is not a float."
        try:
            analyzer.get_unified_regional_analysis(0, 2, 2)
            raise AssertionError("Zero window size was accepted.")
        except ValueError:
            pass
        try:
            analyzer.get_unified_regional_analysis(6, 0, 2)
            raise AssertionError("Zero k-mer length was accepted.")
        except ValueError:
            pass
        try:
            analyzer.get_unified_regional_analysis(6, 2, 0)
            raise AssertionError("Zero repeat length was accepted.")
        except ValueError:
            pass
        print(green("Version 0.1.0 unified regional analysis is online."))
        success += 1
    except Exception as e:
        failure += 1
        print(red(e))
        print(red("Version 0.1.0 unified regional analysis failed."))

    try:
        tests += 1
        ingestor = SequenceIngestor()
        result = ingestor.ingest("atgc\n nnatgc")
        assert isinstance(result, dict), "Sequence ingestion did not return a dictionary."
        assert result["sequence"] == "ATGCNNATGC", "Sequence normalization is incorrect."
        assert result["sequence_type"] == "DNA", "DNA sequence was not detected correctly."
        assert result["length"] == 10, "Ingested sequence length is incorrect."
        assert result["ambiguous_bases"] == 2, "Ambiguous base count is incorrect."
        assert ingestor.get_sequence() == "ATGCNNATGC", "Ingestor did not retain the sequence."
        assert ingestor.get_sequence_type() == "DNA", "Ingestor did not retain the sequence type."
        result = ingestor.ingest("AUGC\nNNNAUGC")
        assert result["sequence"] == "AUGCNNNAUGC", "RNA sequence normalization is incorrect."
        assert result["sequence_type"] == "RNA", "RNA sequence was not detected correctly."
        assert result["length"] == 11, "RNA sequence length is incorrect."
        assert result["ambiguous_bases"] == 3, "RNA ambiguous base count is incorrect."
        result = ingestor.ingest("ACGNNACG")
        assert result["sequence_type"] == "UNKNOWN", "Ambiguous DNA/RNA sequence was not classified as UNKNOWN."
        assert result["length"] == 8, "Unknown sequence length is incorrect."
        try:
            ingestor.ingest("")
            raise AssertionError("Empty sequence was accepted.")
        except ValueError:
            pass
        try:
            ingestor.ingest("ATGCX")
            raise AssertionError("Invalid nucleotide was accepted.")
        except ValueError:
            pass
        try:
            ingestor.ingest(12345)
            raise AssertionError("Non-string sequence was accepted.")
        except TypeError:
            pass
        print(green("Version 0.1.1 raw sequence ingestion is online."))
        success += 1
    except Exception as e:
        failure += 1
        print(red(e))
        print(red("Version 0.1.1 raw sequence ingestion failed."))

    try:
        tests += 1
        from classes.fasta_parser import FastaParser
        parser = FastaParser()
        result = parser.parse(">virus_sample\nATGC\nNNATGC")
        assert isinstance(result, dict), "FASTA parser did not return a dictionary."
        assert result["header"] == "virus_sample", "FASTA header was parsed incorrectly."
        assert result["sequence"] == "ATGCNNATGC", "FASTA sequence was parsed incorrectly."
        assert result["sequence_type"] == "DNA", "FASTA DNA sequence was not detected correctly."
        assert result["length"] == 10, "FASTA sequence length is incorrect."
        assert result["ambiguous_bases"] == 2, "FASTA ambiguous base count is incorrect."
        assert parser.get_header() == "virus_sample", "FASTA parser did not retain the header."
        assert parser.get_sequence() == "ATGCNNATGC", "FASTA parser did not retain the sequence."
        assert parser.get_metadata()["header"] == "virus_sample", "FASTA metadata header is incorrect."
        result = parser.parse(">rna_sample\nAUGC\nNNNAUGC")
        assert result["header"] == "rna_sample", "RNA FASTA header was parsed incorrectly."
        assert result["sequence"] == "AUGCNNNAUGC", "RNA FASTA sequence was parsed incorrectly."
        assert result["sequence_type"] == "RNA", "RNA FASTA sequence was not detected correctly."
        assert result["length"] == 11, "RNA FASTA sequence length is incorrect."
        assert result["ambiguous_bases"] == 3, "RNA FASTA ambiguous base count is incorrect."
        try:
            parser.parse("ATGC")
            raise AssertionError("FASTA content without a header was accepted.")
        except ValueError:
            pass
        try:
            parser.parse(">virus_sample")
            raise AssertionError("FASTA content without a sequence was accepted.")
        except ValueError:
            pass
        try:
            parser.parse(">")
            raise AssertionError("FASTA content with an empty header was accepted.")
        except ValueError:
            pass
        try:
            parser.parse(12345)
            raise AssertionError("Non-string FASTA content was accepted.")
        except TypeError:
            pass
        print(green("Version 0.2.0 FASTA parsing is online."))
        success += 1
    except Exception as e:
        failure += 1
        print(red(e))
        print(red("Version 0.2.0 FASTA parsing failed."))

    try:
        tests += 1
        from classes.configuration import Configuration
        configuration = Configuration()
        assert isinstance(configuration, Configuration), "Configuration did not initialize as a Configuration instance."
        assert configuration.get("window_size") == 100, "Default window size is incorrect."
        assert configuration.get("kmer_size") == 3, "Default k-mer size is incorrect."
        assert configuration.get("repeat_length") == 3, "Default repeat length is incorrect."
        assert configuration.get("motifs") == {}, "Default motif configuration is incorrect."
        configuration.set("window_size", 250)
        assert configuration.get("window_size") == 250, "Updated window size was not stored correctly."
        configuration.set("kmer_size", 5)
        assert configuration.get("kmer_size") == 5, "Updated k-mer size was not stored correctly."
        configuration.set("repeat_length", 4)
        assert configuration.get("repeat_length") == 4, "Updated repeat length was not stored correctly."
        configuration.set("motifs", {"START": "ATG"})
        assert configuration.get("motifs") == {"START": "ATG"}, "Updated motif configuration was not stored correctly."
        settings = configuration.get_all()
        assert isinstance(settings, dict), "Configuration settings did not return a dictionary."
        assert settings["window_size"] == 250, "Configuration dictionary contains an incorrect window size."
        assert settings["kmer_size"] == 5, "Configuration dictionary contains an incorrect k-mer size."
        assert settings["repeat_length"] == 4, "Configuration dictionary contains an incorrect repeat length."
        assert settings["motifs"] == {"START": "ATG"}, "Configuration dictionary contains an incorrect motif configuration."
        settings["window_size"] = 999
        assert configuration.get("window_size") == 250, "Configuration was modified through the returned settings dictionary."
        try:
            configuration.get("unknown")
            raise AssertionError("Unknown configuration setting was accepted.")
        except KeyError:
            pass
        try:
            configuration.set("unknown", 123)
            raise AssertionError("Unknown configuration setting was updated.")
        except KeyError:
            pass
        configuration.reset()
        assert configuration.get("window_size") == 100, "Configuration reset did not restore the default window size."
        assert configuration.get("kmer_size") == 3, "Configuration reset did not restore the default k-mer size."
        assert configuration.get("repeat_length") == 3, "Configuration reset did not restore the default repeat length."
        assert configuration.get("motifs") == {}, "Configuration reset did not restore the default motifs."
        print(green("Version 0.3.0 configuration is online."))
        success += 1
    except Exception as e:
        failure += 1
        print(red(e))
        print(red("Version 0.3.0 configuration failed."))

    try:
        tests += 1
        from classes.feature_normalizer import FeatureNormalizer
        normalizer = FeatureNormalizer()
        value = normalizer.normalize_value(50, 0, 100)
        assert value == 0.5, "Midpoint normalization did not produce 0.5."
        value = normalizer.normalize_value(0, 0, 100)
        assert value == 0.0, "Minimum normalization did not produce 0.0."
        value = normalizer.normalize_value(100, 0, 100)
        assert value == 1.0, "Maximum normalization did not produce 1.0."
        value = normalizer.normalize_value(25, 0, 100, True)
        assert value == 0.75, "Inverted normalization produced an incorrect value."
        value = normalizer.normalize_value(-50, 0, 100)
        assert value == 0.0, "Values below the normalization range were not clamped to 0.0."
        value = normalizer.normalize_value(150, 0, 100)
        assert value == 1.0, "Values above the normalization range were not clamped to 1.0."
        value = normalizer.normalize_value(50, 50, 50)
        assert value == 0.0, "Identical minimum and maximum values did not produce 0.0."
        features = {
            "gc_content": 0.5,
            "complexity": 0.25,
            "repeat_density": 0.75
        }
        ranges = {
            "gc_content": {
                "minimum": 0.0,
                "maximum": 1.0
            },
            "complexity": {
                "minimum": 0.0,
                "maximum": 1.0
            },
            "repeat_density": {
                "minimum": 0.0,
                "maximum": 1.0
            }
        }
        normalized = normalizer.normalize_features(
            features,
            ranges,
            {"complexity"}
        )
        assert isinstance(normalized, dict), "Normalized features did not return a dictionary."
        assert normalized["gc_content"] == 0.5, "GC content normalization is incorrect."
        assert normalized["complexity"] == 0.75, "Inverted complexity normalization is incorrect."
        assert normalized["repeat_density"] == 0.75, "Repeat density normalization is incorrect."
        regional_features = [
            {
                "region": 1,
                "features": {
                    "gc_content": 0.25,
                    "complexity": 0.75
                }
            },
            {
                "region": 2,
                "features": {
                    "gc_content": 0.75,
                    "complexity": 0.25
                }
            }
        ]
        regional_ranges = {
            "gc_content": {
                "minimum": 0.0,
                "maximum": 1.0
            },
            "complexity": {
                "minimum": 0.0,
                "maximum": 1.0
            }
        }
        normalized_regions = normalizer.normalize_regional_features(
            regional_features,
            regional_ranges,
            {"complexity"}
        )
        assert isinstance(normalized_regions, list), "Regional normalization did not return a list."
        assert len(normalized_regions) == 2, "Regional normalization returned an incorrect number of regions."
        assert normalized_regions[0]["region"] == 1, "First normalized region number is incorrect."
        assert normalized_regions[1]["region"] == 2, "Second normalized region number is incorrect."
        assert normalized_regions[0]["features"]["gc_content"] == 0.25, "First region GC normalization is incorrect."
        assert normalized_regions[0]["features"]["complexity"] == 0.25, "First region inverted complexity normalization is incorrect."
        assert normalized_regions[1]["features"]["gc_content"] == 0.75, "Second region GC normalization is incorrect."
        assert normalized_regions[1]["features"]["complexity"] == 0.75, "Second region inverted complexity normalization is incorrect."
        try:
            normalizer.normalize_value("50", 0, 100)
            raise AssertionError("Non-numeric normalization value was accepted.")
        except TypeError:
            pass
        try:
            normalizer.normalize_value(50, 100, 0)
            raise AssertionError("Reversed normalization range was accepted.")
        except ValueError:
            pass
        try:
            normalizer.normalize_features(
                {"gc_content": 0.5},
                {}
            )
            raise AssertionError("Feature without a normalization range was accepted.")
        except KeyError:
            pass
        print(green("Version 0.4.0 feature normalization is online."))
        success += 1
    except Exception as e:
        failure += 1
        print(red(e))
        print(red("Version 0.4.0 feature normalization failed."))

    try:
        tests += 1
        from classes.risk_scorer import RiskScorer
        scorer = RiskScorer({
            "complexity": 1,
            "repeat_density": 2,
            "motif_density": 1
        })
        assert isinstance(scorer, RiskScorer), "Risk scorer did not initialize as a RiskScorer instance."
        assert scorer.get_weights()["complexity"] == 1, "Initial complexity weight is incorrect."
        assert scorer.get_weights()["repeat_density"] == 2, "Initial repeat-density weight is incorrect."
        assert scorer.get_weights()["motif_density"] == 1, "Initial motif-density weight is incorrect."
        normalized_weights = scorer.normalize_weights()
        assert abs(normalized_weights["complexity"] - 0.25) < 0.000001, "Normalized complexity weight is incorrect."
        assert abs(normalized_weights["repeat_density"] - 0.5) < 0.000001, "Normalized repeat-density weight is incorrect."
        assert abs(normalized_weights["motif_density"] - 0.25) < 0.000001, "Normalized motif-density weight is incorrect."
        features = {
            "complexity": 0.8,
            "repeat_density": 0.6,
            "motif_density": 0.4
        }
        result = scorer.score_features(features)
        assert isinstance(result, dict), "Risk scoring did not return a dictionary."
        assert "score" in result, "Risk score is missing from the scoring result."
        assert "contributions" in result, "Feature contributions are missing from the scoring result."
        assert abs(result["score"] - 0.6) < 0.000001, "Combined risk score is incorrect."
        assert abs(result["contributions"]["complexity"] - 0.2) < 0.000001, "Complexity contribution is incorrect."
        assert abs(result["contributions"]["repeat_density"] - 0.3) < 0.000001, "Repeat-density contribution is incorrect."
        assert abs(result["contributions"]["motif_density"] - 0.1) < 0.000001, "Motif-density contribution is incorrect."
        assert abs(sum(result["contributions"].values()) - result["score"]) < 0.000001, "Feature contributions do not sum to the risk score."
        regional_features = [
            {
                "region": 1,
                "features": {
                    "complexity": 0.8,
                    "repeat_density": 0.6,
                    "motif_density": 0.4
                }
            },
            {
                "region": 2,
                "features": {
                    "complexity": 0.2,
                    "repeat_density": 0.4,
                    "motif_density": 0.8
                }
            }
        ]
        scored_regions = scorer.score_regions(regional_features)
        assert isinstance(scored_regions, list), "Regional risk scoring did not return a list."
        assert len(scored_regions) == 2, "Regional risk scoring returned an incorrect number of regions."
        assert scored_regions[0]["region"] == 1, "First scored region number is incorrect."
        assert scored_regions[1]["region"] == 2, "Second scored region number is incorrect."
        assert abs(scored_regions[0]["risk_score"] - 0.6) < 0.000001, "First regional risk score is incorrect."
        assert abs(scored_regions[1]["risk_score"] - 0.45) < 0.000001, "Second regional risk score is incorrect."
        assert "contributions" in scored_regions[0], "First region feature contributions are missing."
        assert "contributions" in scored_regions[1], "Second region feature contributions are missing."
        try:
            scorer.set_weights({})
            raise AssertionError("Empty risk weights were accepted.")
        except ValueError:
            pass
        try:
            scorer.set_weights({"complexity": -1})
            raise AssertionError("Negative risk weight was accepted.")
        except ValueError:
            pass
        try:
            scorer.score_features({"complexity": 0.5})
            raise AssertionError("Incomplete feature set was accepted for risk scoring.")
        except KeyError:
            pass
        try:
            scorer.score_features({
                "complexity": "high",
                "repeat_density": 0.5,
                "motif_density": 0.5
            })
            raise AssertionError("Non-numeric risk feature was accepted.")
        except TypeError:
            pass
        print(green("Version 0.5.0 mutation-risk scoring is online."))
        success += 1
    except Exception as e:
        failure += 1
        print(red(e))
        print(red("Version 0.5.0 mutation-risk scoring failed."))

    print()
    print("===================================")
    print("VIRAL MUTATION ANALYZER TEST SUITE")
    print("===================================")
    print(f"Tests:   {tests}")
    print(green(f"Success: {success}"))
    print(red(f"Failure: {failure}"))
    print("===================================")

    if failure == 0:
        print(green("ALL TESTS PASSED."))
    else:
        print(red("TEST SUITE FAILED."))