import array
from pathlib import Path
import sys
import tempfile
import unittest

from ms2exp.validate.legacy import check_lengths


class LengthValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.truth = self.root/'truth.u32'
        values = array.array('I', [3, 2, 1, 1, 0])
        if sys.byteorder != 'little':
            values.byteswap()
        self.truth.write_bytes(values.tobytes())
        self.output = self.root/'output.lengths'

    def check(self, text, baseline=False):
        self.output.write_text(text)
        return check_lengths(self.output, self.truth, [0, 3, 5], baseline)

    def test_complete_and_baseline(self):
        result = self.check('>0\n3 2\n1\n>q1\n1 0\n')
        self.assertEqual((result['reads'], result['positions'], result['failed']), (2, 5, 0))
        self.assertTrue(result['every_position_checked'])
        self.assertEqual(self.check('>q0\n3 2 1\n', True)['positions'], 3)

    def test_malformed_or_incomplete_records(self):
        cases = ['', '3 2 1\n', '>q0\n3 2\n', '>q0\n3 2 1\n',
                 '>q0\n3 2 1 0\n>q1\n1 0\n', '>q1\n3 2 1\n>q0\n1 0\n',
                 '>q0\n3 2 1\n>q0\n1 0\n', '>q0\n3 1 1\n>q1\n1 0\n',
                 '>q0\n3 -2 1\n>q1\n1 0\n', '>q0\n3 2 2\n>q1\n1 0\n',
                 '>q0\n3 2 1\n>q1\n1 0\n>q2\n0\n']
        for text in cases:
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.check(text)
        with self.assertRaises(ValueError):
            self.check('>q0\n3 2 1\n>q1\n1 0\n', True)

    def test_invalid_truth_dimensions(self):
        self.truth.write_bytes(b'\0'*16)
        with self.assertRaises(ValueError):
            self.check('>q0\n3 2 1\n>q1\n1 0\n')


if __name__ == '__main__':
    unittest.main()
