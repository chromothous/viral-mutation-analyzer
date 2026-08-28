import os
import sys

from classes.analyzer import Analyzer
from classes.logger import Logger

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
        print(red(e))
        print(red("Version 0.0.4 sequence statistics failed."))

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