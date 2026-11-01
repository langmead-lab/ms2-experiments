#!/usr/bin/env python3
"""Stream the source once; samtools selects reads. Never spool the full BAM."""
import datetime, gzip, hashlib, json, pathlib, platform, shutil, subprocess, sys, time
root=pathlib.Path(sys.argv[1]); root.mkdir(parents=True,exist_ok=True)
sam=shutil.which('samtools')
if not sam: raise SystemExit('samtools must be on PATH')
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
m={'status':'running','started_utc':now(),'source_url':'https://ont-open-data.s3.amazonaws.com/gm24385_2023.12/all_pass.vhg002v1.bam','expected_source_bytes':96358011601,'documentation_url':'https://epi2me.nanoporetech.com/gm24385_ncm23_preview/','fraction':0.008,'seed':20261003,'excluded_flags':'0x900 (secondary and supplementary); unmapped retained','samtools_binary':sam,'samtools_sha256':sha(sam),'samtools_version':subprocess.check_output([sam,'--version'],text=True),'python':sys.version,'host':platform.node(),'script_sha256':sha(__file__),'commands':[]}
def save():
    m['updated_utc']=now(); p=root/'manifest.json.tmp'; p.write_text(json.dumps(m,indent=2)+'\n'); p.replace(root/'manifest.json')
def cmd(a): m['commands'].append(a); save(); return a
save()
try:
    headers=subprocess.check_output(cmd(['curl','-fsSI',m['source_url']]),text=True)
    (root/'source-head.txt').write_text(headers)
    etags=[line.split(':',1)[1].strip() for line in headers.splitlines() if line.lower().startswith('etag:')]
    if not etags: raise ValueError('No source ETag')
    m['source_etag']=etags[-1]; m['etag_is_md5']=False
    with (root/'curl.stderr').open('w') as ce, (root/'samtools-view.stderr').open('w') as se:
        fetch=subprocess.Popen(cmd(['curl','--fail','--location','--connect-timeout','30','--speed-limit','1024','--speed-time','300','--header','If-Match: '+etags[-1],'--dump-header',str(root/'source-get-headers.txt'),m['source_url']]),stdout=subprocess.PIPE,stderr=ce)
        select=subprocess.Popen(cmd([sam,'view','--subsample','0.008','--subsample-seed','20261003','-F','0x900','-b','-o',str(root/'sample.bam.partial'),'-']),stdin=subprocess.PIPE,stderr=se)
        h=hashlib.sha256(); n=0; last=time.monotonic()
        try:
            while True:
                b=fetch.stdout.read(8*1024*1024)
                if not b: break
                h.update(b); n+=len(b); select.stdin.write(b)
                if time.monotonic()-last>30:
                    m['source_bytes_streamed']=n; save(); last=time.monotonic()
            select.stdin.close(); fetch.stdout.close()
            fr=fetch.wait(); sr=select.wait()
            if fr or sr: raise RuntimeError(f'curl={fr}, samtools view={sr}')
        finally:
            for p in (fetch,select):
                if p.poll() is None: p.terminate(); p.wait()
    if n != m['expected_source_bytes']: raise ValueError(f'Truncated source: {n}')
    m.update(source_bytes_streamed=n,source_sha256=h.hexdigest()); save()
    bam=root/'sample.bam'; (root/'sample.bam.partial').rename(bam)
    subprocess.run(cmd([sam,'quickcheck','-v',str(bam)]),check=True)
    (root/'sample-header.sam').write_text(subprocess.check_output(cmd([sam,'view','-H',str(bam)]),text=True))
    counts={'primary_records':0,'unmapped':0,'hard_clipped':0,'missing_sequence':0,'missing_quality':0}; names=set()
    p=subprocess.Popen(cmd([sam,'view',str(bam)]),stdout=subprocess.PIPE,text=True)
    for line in p.stdout:
        fields=line.rstrip('\n').split('\t'); name=fields[0]; flag=int(fields[1])
        if name in names or flag & 0x900: raise ValueError('Duplicate or non-primary selected record')
        names.add(name); counts['primary_records']+=1; counts['unmapped']+=bool(flag&4)
        counts['hard_clipped']+=('H' in fields[5]); counts['missing_sequence']+=(fields[9]=='*'); counts['missing_quality']+=(fields[10]=='*')
    if p.wait(): raise RuntimeError('samtools audit failed')
    m['selected_record_audit']=counts; save()
    if counts['hard_clipped'] or counts['missing_sequence'] or counts['missing_quality']: raise ValueError('Selected records cannot faithfully reconstruct complete FASTQ; inspect audit')
    with (root/'sample.fastq.gz.partial').open('wb') as out, (root/'samtools-fastq.stderr').open('w') as err:
        p=subprocess.Popen(cmd([sam,'fastq','-n',str(bam)]),stdout=subprocess.PIPE,stderr=err)
        z=subprocess.Popen(cmd(['gzip','-n','-c']),stdin=p.stdout,stdout=out); p.stdout.close()
        zr=z.wait(); pr=p.wait()
        if zr or pr: raise RuntimeError(f'FASTQ conversion failed {pr}/{zr}')
    total=bases=0; ids=set(); fqhash=hashlib.sha256()
    with gzip.open(root/'sample.fastq.gz.partial','rb') as f, (root/'read-names.txt').open('w') as out:
        while True:
            header=f.readline()
            if not header: break
            seq=f.readline(); plus=f.readline(); qual=f.readline()
            for line in (header,seq,plus,qual): fqhash.update(line)
            if not header.startswith(b'@') or not plus.startswith(b'+') or len(seq.rstrip())!=len(qual.rstrip()): raise ValueError('Invalid FASTQ')
            name=header[1:].split()[0].decode()
            if name in ids: raise ValueError('Duplicate FASTQ name')
            ids.add(name); out.write(name+'\n'); total+=1; bases+=len(seq.rstrip())
    if ids != names: raise ValueError('FASTQ membership mismatch')
    (root/'sample.fastq.gz.partial').rename(root/'sample.fastq.gz')
    m.update(status='complete',reads=total,bases=bases,uncompressed_fastq_sha256=fqhash.hexdigest(),gzip_version=subprocess.check_output(['gzip','--version'],text=True),files=[{'name':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in [bam,root/'sample.fastq.gz',root/'read-names.txt']])
except Exception as e:
    m.update(status='failed',error=repr(e)); raise
finally: save()
