"""Exercise variable lengths, output reordering, missing reads and mismatches."""
import struct

import pytest

from ms2exp.validate.ms import check


def test_variable_lengths_reordered_outputs_and_malformed_outputs(tmp_path):
 d=tmp_path;truth=d/'truth';offsets=[0,3,5];row={'stdout':str(d/'rope')}
 (d/'rope').write_text('q0\t0\t3\t1\nq1\t0\t2\t1\n')
 assert check('ropebwt3',d,row,truth,offsets,True)['positions']==5
 (d/'output').write_text('q1\t2 1\nq0\t3 2 1\n')
 check('Orbit-ms',d,row,truth,offsets)
 with (d/'output.len').open('wb') as f:
  for n,values in [('q1',[2,1]),('q0',[3,2,1])]:f.write(('>'+n+'\n').encode()+struct.pack('<Q',len(values))+struct.pack('<'+'I'*len(values),*values)+b'\n')
 check('TeraMS-dual',d,row,truth,offsets)
 for bad in ['q0\t3 2 1\n','q0\t3 1 1\nq1\t2 1\n','q0\t3 2 1\nq0\t3 2 1\n']:
  (d/'output').write_text(bad)
  with pytest.raises((ValueError,AssertionError)):check('Orbit-ms',d,row,truth,offsets)
