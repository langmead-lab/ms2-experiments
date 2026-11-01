"""Shared pieces of the ms2 freeze benchmark (BENCH_FREEZE_DESIGN.md): manifest access, hashing and stores.

A store is either a local directory or an S3 prefix (s3://bucket/prefix). Inputs live in an
inputs store under their logical names; results go to a results store. The same code runs on
Rockfish (local stores on /vast or /scratch) and on AWS (S3 stores, through the aws CLI so the
boxes need no Python packages beyond the standard library).
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

MANIFEST_FORMAT = 'ms2-freeze-manifest-v1'
RESULT_FORMAT = 'ms2-freeze-cell-v1'
# The repository root: the code directory a manifest's build steps name as {code}.
REPO = Path(__file__).resolve().parents[1]


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(16 << 20), b''):
            h.update(block)
    return h.hexdigest()


def json_sha256(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def load_manifest(path):
    manifest = json.loads(Path(path).read_text())
    if manifest.get('format') != MANIFEST_FORMAT:
        raise ValueError('Not a freeze manifest: ' + str(path))
    return manifest


def find_job(manifest, job_id):
    jobs = [j for j in manifest['jobs'] if j['id'] == job_id]
    if len(jobs) != 1:
        raise KeyError('No single job ' + job_id)
    return jobs[0]


def cell_key(cell_id):
    """Directory name of a cell id (group/fig_id)."""
    return cell_id.replace('/', '__')


class Store:
    """A directory or an S3 prefix holding files under relative names."""

    def __init__(self, root):
        self.root = str(root).rstrip('/')
        self.s3 = self.root.startswith('s3://')
        if not self.s3:
            self.path = Path(self.root)

    def __repr__(self):
        return f'Store({self.root!r})'

    def url(self, name):
        return f'{self.root}/{name}'

    def _aws(self, *args, check=True, capture=True):
        argv = ['aws', *args]
        if os.environ.get('FREEZE_S3_NO_SIGN'):
            argv.append('--no-sign-request')
        result = subprocess.run(argv, capture_output=capture, text=True)
        if check and result.returncode:
            raise RuntimeError(f'{" ".join(argv)}: exit {result.returncode}: {result.stderr.strip()[-2000:]}')
        return result

    def exists(self, name):
        if not self.s3:
            return (self.path/name).exists()
        bucket, key = self.url(name)[5:].split('/', 1)
        return self._aws('s3api', 'head-object', '--bucket', bucket, '--key', key, check=False).returncode == 0

    def get(self, name, dest):
        """Copy one file to dest (a local path)."""
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not self.s3:
            shutil.copyfile(self.path/name, dest)
        else:
            self._aws('s3', 'cp', '--only-show-errors', self.url(name), str(dest))

    def local_path(self, name):
        """The file itself when the store is a directory; None for S3."""
        return None if self.s3 else self.path/name

    def put(self, src, name):
        if not self.s3:
            dest = self.path/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            tmp = dest.with_name(dest.name + '.tmp')
            shutil.copyfile(src, tmp)
            tmp.replace(dest)
        else:
            self._aws('s3', 'cp', '--only-show-errors', str(src), self.url(name))

    def put_tree(self, src, name):
        src = Path(src)
        for path in sorted(p for p in src.rglob('*') if p.is_file()):
            self.put(path, f'{name}/{path.relative_to(src)}')

    def put_json(self, value, name):
        with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
            f.write(json.dumps(value, indent=2) + '\n')
        try:
            self.put(f.name, name)
        finally:
            os.unlink(f.name)

    def get_json(self, name):
        if not self.s3:
            return json.loads((self.path/name).read_text())
        with tempfile.TemporaryDirectory() as d:
            self.get(name, Path(d)/'x.json')
            return json.loads((Path(d)/'x.json').read_text())

    def put_if_absent(self, value, name):
        """Create name holding value unless it exists; True when this call created it."""
        if not self.s3:
            dest = self.path/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            try:
                with dest.open('x') as f:
                    f.write(json.dumps(value, indent=2) + '\n')
                return True
            except FileExistsError:
                return False
        bucket, key = self.url(name)[5:].split('/', 1)
        with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as f:
            f.write(json.dumps(value, indent=2) + '\n')
        try:
            # S3 conditional writes make the claim atomic across boxes.
            result = self._aws('s3api', 'put-object', '--bucket', bucket, '--key', key, '--body', f.name,
                               '--if-none-match', '*', check=False)
        finally:
            os.unlink(f.name)
        if result.returncode == 0:
            return True
        if 'PreconditionFailed' in result.stderr or 'ConditionalRequestConflict' in result.stderr:
            return False
        raise RuntimeError('Conditional put failed: ' + result.stderr.strip()[-2000:])

    def names(self, prefix):
        """Every file name under prefix (relative to the store root), in one listing."""
        prefix = prefix.rstrip('/')
        if not self.s3:
            base = self.path/prefix
            return [str(p.relative_to(self.path)) for p in base.rglob('*') if p.is_file()] if base.exists() else []
        result = self._aws('s3', 'ls', '--recursive', self.url(prefix) + '/', check=False)
        _, key_root = (self.root[5:] + '/').split('/', 1)
        return [line.split(None, 3)[3][len(key_root):] for line in result.stdout.splitlines()
                if len(line.split(None, 3)) == 4]

    def delete(self, name):
        if not self.s3:
            (self.path/name).unlink(missing_ok=True)
        else:
            self._aws('s3', 'rm', '--only-show-errors', self.url(name), check=False)


def now():
    return time.time()
