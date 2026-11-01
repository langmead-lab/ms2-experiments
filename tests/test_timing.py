"""Phase timers of ms2 and ropebwt3 under the load-query-separate policy, and counting-build stats lines."""
import unittest

from ms2exp.run.timing import orbit, ropebwt3, parse, parse_stats

ORBIT = ('smem-all: positions as contig:strand:offset, 2 sequences\n'
         'ftab: k=14 entry_bits=47 bytes=10 load_s=0.7\n'
         'batch: reads=2 bases=4 threads=1 index_load_s=59.3 ftab_k=14 ftab_s=0.7 query_s=3.7 '
         'read_s=0.01 parse_s=0.05 format_s=0.9 orient=fwd flipped=0 orient_s=0 total_s=4.75\n')
ROPE = ('[M::rb3_fmi_load_all::6.796*1.00] loaded the BWT\n'
        '[M::rb3_fmi_load_all::11.628*1.00] loaded the sequence names and lengths\n'
        '[M::worker_pipeline::50.000*1.00] processed 3 sequences\n'
        '[M::worker_pipeline::61.628*1.00] processed 2 sequences\n'
        '[M::main] Real time: 61.700 sec; CPU: 61.600 sec; Peak RSS: 24.7 GB\n')
QUERY_LINE = ('query: reads=2 bases=9 threads=1 index_load_s=210.5 ftab_k=14 ftab_s=1.5 query_s=2.25 '
              'read_s=0.01 parse_s=0.02 format_s=0.3 orient=fwd flipped=0 orient_s=0 total_s=3.5 '
              'query_ns_per_base=250 wall_ns_per_base=388.9\n')
COUNT_LINES = ('count: backward search table 1932735283 bytes, row shift 6, built in 17.6 s\n'
               'count: method=adaptive budget=8 budget_per_len=0.25 short_len=0 smems=162700 walked=1000 '
               'switched=20 searched=161680 walk_steps=5000 search_steps=900000 top_moves=1 bot_moves=2 '
               'scan_rows=3 lf_ff=4 far_ends=5 sum_rows=6 list_steps=7 walk_s=0.01 search_s=1.2 list_s=0.4 '
               'count_init_s=17.6\n')
TMS_STATS = ('stats: bases=35 repositions/base=0.1 psi_steps/base=1.5 scan_rows/rep=2 len/rep=3 phi_steps/base=0.25 '
             'dist/rep=1 dist1_frac=0.5 lce/rep=2 capped_frac=0 scan_visits/rep=1 lf_ff/step=0.2 walk_visits/base=0.01\n'
             'stats: bases=35 repositions/base=0 lf_ff/step=-nan walk_up/rep=-nan walk_down/rep=-nan\n')
MS_STATS = ('stats: bases=0 repositions/base=-nan psi_steps/base=-nan scan_rows/rep=-nan len/rep=-nan phi_steps/base=-nan '
            'dist/rep=-nan dist1_frac=-nan lce/rep=-nan capped_frac=-nan scan_visits/rep=-nan lf_ff/step=-nan walk_visits/base=-nan\n'
            'stats: bases=35 repositions/base=0.2 lf_ff/step=0.1 walk_up/rep=3 walk_down/rep=2\n')
DUALSUM_STATS = (TMS_STATS + 'stats: dualsum up_raced=4 down_raced=2 up_psi_steps=5 down_psi_steps=6 psi_decided=3 '
                 'phi_decided=1 psi_steps/base=1.5\n'
                 'stats: dualsum up_settled=7 up_walked=1 down_settled=2 down_walked=0 partial_steps=9 cont_walks=1 '
                 'cont_steps=2 partial_steps/rep=0.5 cont_steps/rep=0.25\n')


class ParserTests(unittest.TestCase):
    def test_orbit_load_is_index_plus_ftab_and_query_is_query_plus_format(self):
        phase = orbit(ORBIT)
        self.assertAlmostEqual(phase['load_seconds'], 60.0)
        self.assertAlmostEqual(phase['query_seconds'], 4.75)
        self.assertEqual(phase['components']['total_s'], 4.75)
        with self.assertRaises(ValueError):
            orbit(ORBIT.replace('load_s=0.7', 'load_s=0.9'))
        with self.assertRaises(ValueError):
            orbit('smem-all: no timer\n')

    def test_ropebwt3_uses_last_loaded_and_last_processed_stamps(self):
        phase = ropebwt3(ROPE)
        self.assertAlmostEqual(phase['load_seconds'], 11.628)
        self.assertAlmostEqual(phase['query_seconds'], 50.0)
        self.assertEqual(phase['components']['processed_sequences'], 5)
        with self.assertRaises(ValueError):
            ropebwt3(ROPE.replace('processed', 'skipped'))

    def test_parse_dispatches_by_tool(self):
        self.assertEqual(parse('orbit', ORBIT), orbit(ORBIT))
        self.assertEqual(parse('ropebwt3', ROPE), ropebwt3(ROPE))


class QueryLineTests(unittest.TestCase):
    def test_query_summary_line(self):
        phase = orbit('smem-all: positions as contig:strand:offset, 2 sequences\n' + COUNT_LINES + QUERY_LINE)
        self.assertEqual((phase['load_seconds'], phase['query_seconds']), (212.0, 3.5))
        self.assertEqual(phase['source'], 'orbit-query-line')

    def test_query_error_lines_are_not_summaries(self):
        phase = orbit('query: needs a tms index\n' + QUERY_LINE)
        self.assertEqual(phase['query_seconds'], 3.5)
        with self.assertRaises(ValueError):
            orbit('query: needs a tms index\n')
        with self.assertRaises(ValueError):
            orbit(QUERY_LINE + QUERY_LINE)


class StatsTests(unittest.TestCase):
    def test_parse(self):
        tms = parse_stats(TMS_STATS, 'psi')['ops_per_base']
        self.assertEqual((tms['psi'], tms['phi'], tms['total'], tms['repositions']), (1.5, 0.25, 2.75, 0.1))
        self.assertAlmostEqual(tms['scan_rows_per_base'], 0.2)
        ms = parse_stats(MS_STATS, 'ms')['ops_per_base']
        self.assertAlmostEqual(ms['walk'], 1.0)
        self.assertAlmostEqual(ms['total'], 2.0)
        self.assertEqual((ms['psi'], ms['phi']), (0.0, 0.0))
        with self.assertRaises(ValueError):
            parse_stats('stats: bases=3\n', 'dual')

    def test_summary_lines(self):
        stats = parse_stats(DUALSUM_STATS, 'dualsum')
        self.assertEqual(stats['ops_per_base']['total'], 2.75)
        self.assertEqual(stats['summary_lines']['dualsum']['psi_decided'], 3)
        self.assertEqual(stats['summary_lines']['dualsum']['partial_steps/rep'], 0.5)
        self.assertNotIn('summary_lines', parse_stats(TMS_STATS, 'psi'))


if __name__ == '__main__':
    unittest.main()
