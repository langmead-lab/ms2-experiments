#!/usr/bin/env python3
"""Resumable ENA source downloads; one transfer, atomic publication after MD5."""
import concurrent.futures,datetime,fcntl,hashlib,json,os,subprocess,sys,time
from pathlib import Path
root=Path(sys.argv[1]).resolve()
(root/'logs').mkdir(parents=True,exist_ok=True)
lock=(root/'logs/download.lock').open('w')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
files=json.loads((root/'manifests/ena-files.json').read_text())['files']
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(path,data):
 temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(path)
def fetch(f):
 dst=root/f['path'];dst.parent.mkdir(parents=True,exist_ok=True)
 part=dst.with_name(dst.name+'.part');status=root/'logs'/(dst.name+'.status.json')
 if dst.exists() and dst.stat().st_size==f['bytes'] and status.exists():
  prior=json.loads(status.read_text())
  if prior.get('state')=='verified' and prior.get('observed_md5')==f['md5']:
   print(now(),'ALREADY VERIFIED',f['path'],flush=True);return True
 state=dict(f,started_utc=now(),state='downloading')
 def update(**kw):state.update(kw);write(status,state)
 try:
  if not dst.exists():
   with (root/'logs'/(dst.name+'.curl.log')).open('a') as log:
    cmd=['curl','--fail','--location','--continue-at','-','--retry','20','--retry-connrefused','--retry-delay','30','--connect-timeout','30','--speed-limit','1024','--speed-time','180','--output',str(part),f['url']]
    update(command=cmd)
    attempt=0;stalls=0;deadline=time.monotonic()+48*3600
    while True:
     attempt+=1
     before=part.stat().st_size if part.exists() else 0
     if before==f['bytes']:break
     if time.monotonic()>deadline:raise TimeoutError('48-hour per-file transfer budget exhausted')
     update(state='downloading',attempt=attempt)
     result=subprocess.run(cmd,stdout=log,stderr=log)
     after=part.stat().st_size if part.exists() else 0
     if result.returncode==0:break
     progressed=after>before
     stalls=0 if progressed else stalls+1
     if stalls>=6:raise subprocess.CalledProcessError(result.returncode,cmd)
     delay=15 if progressed else min(60*2**(stalls-1),900)
     update(state='retry_wait',last_curl_exit=result.returncode,retry_delay_seconds=delay,consecutive_stalls=stalls,downloaded_bytes=after)
     print(now(),'RETRY',f['path'],'curl exit',result.returncode,'new bytes',after-before,'wait',delay,flush=True)
     time.sleep(delay)
   candidate=part
  else:candidate=dst
  update(state='verifying',download_finished_utc=now())
  if candidate.stat().st_size!=f['bytes']:raise ValueError('File size does not match ENA')
  h=hashlib.md5()
  with candidate.open('rb') as inp:
   for chunk in iter(lambda:inp.read(8*1024*1024),b''):h.update(chunk)
  if h.hexdigest()!=f['md5']:raise ValueError('MD5 does not match ENA; file retained for inspection')
  if candidate==part:part.replace(dst)
  update(state='verified',finished_utc=now(),observed_md5=h.hexdigest())
  print(now(),'VERIFIED',f['path'],flush=True);return True
 except Exception as exc:
  update(state='failed',finished_utc=now(),error=str(exc));print(now(),'FAILED',f['path'],str(exc),flush=True);return False
oldexit=root/'logs/download.exit'
if oldexit.exists():oldexit.rename(root/'logs'/('download.exit.previous.'+str(int(time.time()))))
write(root/'logs/campaign.json',dict(started_utc=now(),pid=os.getpid(),files=len(files),expected_bytes=sum(f['bytes'] for f in files)))
def remaining(f):
 p=root/f['path'];q=p.with_name(p.name+'.part')
 return max(0,f['bytes']-(p.stat().st_size if p.exists() else q.stat().st_size if q.exists() else 0))
files.sort(key=remaining)
with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:ok=list(pool.map(fetch,files))
code=0 if all(ok) else 1
(root/'logs/download.exit').write_text(str(code)+'\n')
print(now(),'COMPLETE',sum(ok),'of',len(ok),flush=True)
sys.exit(code)
