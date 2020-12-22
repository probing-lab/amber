import glob
import unittest

from mora.input import InputParser
from mora.utils import set_log_level as set_mora_log_level, LOG_NOTHING as MORA_LOG_NOTHING
from src import decide_termination
from src.utils import Answer, set_log_level, LOG_NOTHING

benchmarks_past = glob.glob("tests/benchmarks/past/*")
benchmarks_ast = glob.glob("tests/benchmarks/ast/*")
benchmarks_nast = glob.glob("tests/benchmarks/non-ast/*")

def test_generator_past(benchmark):
    def test(self):
        input_parser = InputParser()
        input_parser.set_source(benchmark)
        program = input_parser.parse_source()
        result = decide_termination(program)
        self.assertEqual(result.PAST, Answer.TRUE, benchmark)
        self.assertEqual(result.AST, Answer.TRUE, benchmark)
        self.assertTrue(len(result.witnesses) > 0, benchmark)
    return test


def test_generator_ast(benchmark):
    def test(self):
        input_parser = InputParser()
        input_parser.set_source(benchmark)
        program = input_parser.parse_source()
        result = decide_termination(program)
        self.assertEqual(result.AST, Answer.TRUE, benchmark)
        self.assertTrue(len(result.witnesses) > 0, benchmark)
    return test


def test_generator_nast(benchmark):
    def test(self):
        input_parser = InputParser()
        input_parser.set_source(benchmark)
        program = input_parser.parse_source()
        result = decide_termination(program)
        self.assertEqual(result.PAST, Answer.FALSE, benchmark)
        self.assertEqual(result.AST, Answer.FALSE, benchmark)
        self.assertTrue(len(result.witnesses) > 0, benchmark)
    return test


class TestBenchmarks(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        set_mora_log_level(MORA_LOG_NOTHING)
        set_log_level(LOG_NOTHING)


for b in benchmarks_past:
    setattr(TestBenchmarks, f"test_past_{b}", test_generator_past(b))

for b in benchmarks_ast:
    setattr(TestBenchmarks, f"test_ast_{b}", test_generator_ast(b))

for b in benchmarks_nast:
    setattr(TestBenchmarks, f"test_nast_{b}", test_generator_nast(b))


if __name__ == '__main__':
    unittest.main()
