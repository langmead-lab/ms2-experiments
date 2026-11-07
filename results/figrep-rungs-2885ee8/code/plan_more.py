#!/usr/bin/env python3
"""Larger chr22 cells (Fig Rep series B above 100 haplotypes) by plan.py's nested-prefix rule.

plan.py froze a panel of the first 100 haplotypes of r2_assemblies_only.txt (CHM13 and GRCh38
removed). This extends the panel by the same rule to every haplotype in that list, maps chr22
from each haplotype's RagTag AGP exactly as plan.py does, checks that the first 100 entries
reproduce the frozen panel, and writes chr22 cells for the counts given to OUT.
"""
import argparse,hashlib,json,re
from pathlib import Path
BASE=Path('/data/blangme2/fasta/hprc-v2')
DICT=Path('/data/blangme2/movi/hprc-yr2/index/hprc-yr2.seqdict.json')
DATA=Path('/scratch16/blangme2/ms-bakery/data/ladder')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('out');ap.add_argument('counts',nargs='+')
 args=ap.parse_args()
 frozen=json.loads((DATA/'plan.json').read_text())
 seq={h['name']:h for h in json.loads(DICT.read_text())['haplotypes']}
 names=[Path(x).stem for x in (BASE/'r2_assemblies_only.txt').read_text().splitlines()]
 names=[n for n in names if not any(x in n for x in ['CHM13','GRCh38'])]
 assert len(names)==len(set(names))
 panel=[]
 for name in names:
  h=seq[name];contigs={c['name']:c for c in h['contigs']}
  agp=BASE/'ragtag'/h['pansn_prefix'].replace('#','.')/'ragtag_output/ragtag.scaffold.agp'
  mapped=[]
  for line in agp.read_text().splitlines():
   if not line or line.startswith('#'):continue
   f=line.split('\t');match=re.search(r'#(chr(?:[0-9]+|X|Y|M))_RagTag$',f[0])
   if f[4]!='W' or not match or match[1]!='chr22':continue
   c=contigs[f[5]];begin,end=int(f[6])-1,int(f[7])
   assert 0<=begin<end<=c['length']
   mapped.append({'contig':f[5],'begin':begin,'end':end,'orientation':f[8],'contig_md5':c['md5']})
  if not mapped:raise ValueError(f'{name}: no chr22 mapping')
  panel.append({'name':name,'mapping':str(agp),'mapping_sha256':sha(agp),'chr22':mapped})
 for a,b in zip(frozen['panel'],panel):
  assert a['name']==b['name'] and a['chromosomes']['chr22']==b['chr22'],a['name']
 cells=[]
 for k in args.counts:
  count=len(panel) if k=='all' else int(k)
  records=[{'haplotype':h['name'],**r} for h in panel[:count] for r in h['chr22']]
  assert len({(r['haplotype'],r['contig'],r['begin'],r['end']) for r in records})==len(records)
  n=sum(r['end']-r['begin'] for r in records)
  cells.append({'id':f'chr22-h{count}','coverage':'chr22','chromosomes':['chr22'],'haplotypes':count,
                'panel':[h['name'] for h in panel[:count]],'native_bases_before_gap_splitting':n,
                'components':records,'planning_disk_budget_bytes':50*n})
 extra={k:v for k,v in frozen.items() if k not in ('cells','panel')}
 extra.update(cells=cells,derived_from=str(DATA/'plan.json'),derived_from_sha256=sha(DATA/'plan.json'),
              panel_size=len(panel),
              note='Fig Rep series B above 100 haplotypes: the frozen panel extended by plan.py\'s rule to every r2_assemblies_only haplotype (CHM13 and GRCh38 removed), chr22 only')
 out=Path(args.out)
 if out.exists():raise FileExistsError(out)
 out.write_text(json.dumps(extra,indent=2)+'\n')
 for c in cells:print(c['id'],c['native_bases_before_gap_splitting'],len(c['components']),flush=True)
if __name__=='__main__':main()
