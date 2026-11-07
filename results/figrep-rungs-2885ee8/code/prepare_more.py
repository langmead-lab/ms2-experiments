#!/usr/bin/env python3
"""bench/ladder/prepare.py for the plan-more.json cells, forward reference only.

Components are extracted by the same agc getctg command and split at non-ACGT runs with the
same record names as prepare.py. A component already in the shared ladder cache on scratch16
is read from there; new ones go to a cache next to the plan. reference-both.* is not written:
these cells are built by ropebwt3 from reference.fa, as chr22-h10 and chr22-h100 were.
"""
import argparse,fcntl,hashlib,json,re,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
SHARED=Path('/scratch16/blangme2/ms-bakery/data/ladder/components')
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
def records(path):
 name=None;parts=[]
 with path.open() as f:
  for line in f:
   if line.startswith('>'):
    if name is not None:yield name,''.join(parts)
    name=line[1:].strip();parts=[]
   else:parts.append(line.strip())
  if name is not None:yield name,''.join(parts)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('plan');ap.add_argument('cell');ap.add_argument('--jobs',type=int,default=8)
 args=ap.parse_args()
 plan_path=Path(args.plan).resolve();plan=json.loads(plan_path.read_text())
 cell=next(c for c in plan['cells'] if c['id']==args.cell)
 data=plan_path.parent/'ladder';out=data/args.cell;out.mkdir(parents=True,exist_ok=True)
 if (out/'manifest.json').exists():print('Already prepared',args.cell);return
 cache=data/'components';cache.mkdir(exist_ok=True)
 def fetch(r):
  key=hashlib.sha256(json.dumps(r,sort_keys=True).encode()).hexdigest()
  if (SHARED/(key+'.fa')).exists():return SHARED/(key+'.fa')
  path=cache/(key+'.fa')
  with (cache/(key+'.lock')).open('w') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX)
   if not path.exists():
    tmp=path.with_suffix('.tmp')
    cmd=['/home/blangme2/bin/agc','getctg','-t','2',plan['archive'],f"{r['contig']}@{r['haplotype']}:{r['begin']}-{r['end']-1}"]
    with tmp.open('w') as f:subprocess.run(cmd,stdout=f,check=True)
    seq=''.join(s for _,s in records(tmp)).upper()
    if len(seq)!=r['end']-r['begin']:raise ValueError('AGC length mismatch: '+str(r))
    tmp.replace(path)
  return path
 with ThreadPoolExecutor(args.jobs) as pool:paths=list(pool.map(fetch,cell['components']))
 print('components ready',len(paths),flush=True)
 metadata=[];total=0;forward=out/'reference.fa'
 with forward.open('w') as dst:
  for i,(r,path) in enumerate(zip(cell['components'],paths)):
   seq=''.join(s for _,s in records(path)).upper();segments=[]
   for j,m in enumerate(re.finditer('[ACGT]+',seq)):
    name=f'h{r["haplotype"].split("_")[0]}_c{i}_s{j}'
    dst.write(f'>{name}\n{m.group()}\n');total+=len(m.group())
    segments.append({'record':name,'begin':r['begin']+m.start(),'end':r['begin']+m.end()})
   metadata.append({**r,'extracted_fasta_sha256':digest(path),'segments':segments})
 manifest={k:v for k,v in cell.items() if k!='components'}
 manifest.update(snapshot=plan['snapshot'],archive=plan['archive'],plan_sha256=digest(plan_path),plan_file=str(plan_path),
                 sources=metadata,forward_acgt_bases=total,reference_sha256=digest(forward),
                 both_strands=True,boundary_policy=plan['boundary_policy'])
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print('PREPARED',args.cell,total,'forward ACGT bases',flush=True)
if __name__=='__main__':main()
