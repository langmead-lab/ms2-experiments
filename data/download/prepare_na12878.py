#!/usr/bin/env python3
"""Acquire one complete published flowcell; no read sampling or filtering."""
import datetime, gzip, hashlib, json, pathlib, platform, shutil, subprocess, sys, tarfile
root = pathlib.Path(sys.argv[1]); root.mkdir(parents=True, exist_ok=True)
url = 'https://s3.amazonaws.com/nanopore-human-wgs/rel6/FASTQTars/FAB42260-4177064552_Multi.tar'
archive = root / 'FAB42260-4177064552_Multi.tar'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
def run(args):
    print(json.dumps(args), flush=True); subprocess.run(args, check=True)
manifest={'source_url':url, 'source_documentation':'https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md', 'expected_source_bytes':1795297280, 'selection':'All FASTQ records from complete rel6 flowcell FAB42260; no filtering or sampling', 'python':sys.version, 'host':platform.node(), 'script_sha256':sha(pathlib.Path(__file__)), 'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
try:
    partial=archive.with_suffix('.tar.partial')
    if not archive.exists():
        run(['curl','--fail','--location','--retry','8','--retry-delay','15','--connect-timeout','30','--speed-limit','1024','--speed-time','180','--continue-at','-','--dump-header',str(root/'source-http-headers.txt'),'--output',str(partial),url])
        if partial.stat().st_size != manifest['expected_source_bytes']: raise ValueError('source size mismatch')
        partial.rename(archive)
    if archive.stat().st_size != manifest['expected_source_bytes']: raise ValueError('source size mismatch')
    manifest['archive_sha256']=sha(archive)
    manifest['curl_version']=subprocess.check_output(['curl','--version'],text=True)
    out=root/'reads.fastq.gz.partial'; ids=root/'read-names.txt.partial'
    files=[]; count=bases=0; seen=set()
    with tarfile.open(archive) as tar, out.open('wb') as dest, ids.open('w') as names:
        members=sorted((m for m in tar.getmembers() if m.isfile() and m.name.lower().endswith(('.fastq.gz','.fq.gz'))),key=lambda m:m.name)
        if not members: raise ValueError('No gzipped FASTQ members found')
        (root/'tar-members.json').write_text(json.dumps([{'name':m.name,'bytes':m.size} for m in tar.getmembers()],indent=2)+'\n')
        for m in members:
            h=hashlib.sha256()
            with tar.extractfile(m) as source:
                for b in iter(lambda:source.read(8*1024*1024),b''): h.update(b); dest.write(b)
            files.append({'member':m.name,'bytes':m.size,'sha256':h.hexdigest()})
            with tar.extractfile(m) as source, gzip.open(source,'rt') as fq:
                while True:
                    header=fq.readline()
                    if not header: break
                    seq=fq.readline().rstrip('\r\n'); plus=fq.readline(); qual=fq.readline().rstrip('\r\n')
                    if not header.startswith('@') or not plus.startswith('+') or len(seq)!=len(qual): raise ValueError('invalid FASTQ')
                    name=header[1:].split()[0]
                    if name in seen: raise ValueError('duplicate read name: '+name)
                    seen.add(name); names.write(name+'\n'); count+=1; bases+=len(seq)
    out.rename(root/'reads.fastq.gz'); ids.rename(root/'read-names.txt')
    manifest.update(status='complete', reads=count, bases=bases, member_order='lexicographic pathname; concatenate original gzip members without recompression', members=files, files=[{'name':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in [root/'reads.fastq.gz',root/'read-names.txt']])
except Exception as e:
    manifest.update(status='failed',error=repr(e)); raise
finally:
    manifest['updated_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
