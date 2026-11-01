"""Whole-output matching-statistics check: a tool's MS, rebuilt from its output, against a u32 truth file.

records() reads ropebwt3 mem SMEMs (rebuilding MS from them), ms2's tab-separated MS output, or the
binary .len output of other tools, and check() compares every read with the truth file.
"""
import array
import mmap
from pathlib import Path
import struct
import sys

def records(mode,folder,row,offsets):
 if mode=='ropebwt3':
  last=None;values=None
  with Path(row['stdout']).open() as f:
   for line in f:
    if not line.strip():continue
    name,a,b,count,*rest=line.split();i=int(name[1:]);assert name==f'q{i}' and 0<=i<len(offsets)-1
    if name!=last:
     if last is not None:
      for j in range(1,len(values)):values[j]=max(values[j],values[j-1]-1)
      yield last,values
     last=name;values=[0]*(offsets[i+1]-offsets[i])
    a,b,count=int(a),int(b),int(count)
    if count==0:continue
    assert 0<=a<b<=len(values)
    values[a]=max(values[a],b-a)
  if last is not None:
   for j in range(1,len(values)):values[j]=max(values[j],values[j-1]-1)
   yield last,values
 elif mode.startswith('Orbit'):
  with (folder/'output').open() as f:
   for line in f:
    fields=line.rstrip('\n').split('\t');yield fields[0],list(map(int,fields[1].split()))
 else:
  with (folder/'output.len').open('rb') as f:
   while h:=f.readline():
    assert h.startswith(b'>');n=struct.unpack('<Q',f.read(8))[0];a=array.array('I');a.frombytes(f.read(4*n))
    if sys.byteorder!='little':a.byteswap()
    assert len(a)==n and f.read(1)==b'\n';yield h[1:].strip().decode(),a


def check(mode,folder,row,truth,offsets,write=False,one=False):
 count=1 if one else len(offsets)-1;seen=bytearray(count)
 with truth.open('w+b' if write else 'rb') as f:
  if write:f.truncate(offsets[-1]*4)
  with mmap.mmap(f.fileno(),0,access=mmap.ACCESS_WRITE if write else mmap.ACCESS_READ) as mm:
   for name,values in records(mode,folder,row,offsets):
    i=int(name[1:]);assert name==f'q{i}' and 0<=i<count and not seen[i],('duplicate/unexpected read',name)
    assert len(values)==offsets[i+1]-offsets[i],('length',name)
    seen[i]=1;a=array.array('I',values)
    if sys.byteorder!='little':a.byteswap()
    data=a.tobytes();start=offsets[i]*4;end=offsets[i+1]*4
    if write:mm[start:end]=data
    elif mm[start:end]!=data:raise ValueError('MS mismatch '+name)
 assert all(seen),('missing reads',sum(seen),count)
 return {'reads':count,'positions':offsets[count],'failed':0}
