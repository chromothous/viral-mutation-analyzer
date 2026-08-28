import os
import sys

from classes.analyzer import Analyzer

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