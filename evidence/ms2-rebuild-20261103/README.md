# Rebuild of the freeze ms2 binaries, devlangmead1, 2026-11-03

The freeze build job (`build-tools` of `freeze-2885ee8-3c8c8967`) built `ms2-bench`, `ms2-stats`,
`ms2-replog` and `ms_check` on an AWS box under Amazon Linux 2023 on 2026-10-31, and its own build
record stayed on that box. This rebuild ran the same command in an Amazon Linux 2023 container at the
same paths, and produced the four pinned binaries byte for byte (`result.txt`). Its build record is
`results/freeze-2885ee8-3c8c8967/build-tools/build-2885ee8-ms2.rebuilt.json`.

How it ran (Rockfish's singularity 3.11.4, without root):

1. `singularity build --sandbox al2023 docker://amazonlinux:2023` (image
   `al2023-container-2023.12.20261102.0`, config sha256 846d3310af9eaf26b799cc341729ad5f8c81b3fa4e05015ef88cf3edb5fec5a2; its packages are in
   `image-rpms.txt`).
2. `download_rpms.py gcc-c++ libstdc++-static zlib-devel make git python3`, run in the container,
   resolved and downloaded the build packages the boxes install (`dnf install git make gcc-c++ libstdc++-static
   zlib-devel`, without the boxes' python3-numpy), which were then unpacked into the sandbox with `rpm2cpio | cpio`
   (`added-rpms.txt`): GCC 11.5.0-5.amzn2023.0.5, glibc-devel 2.34-231.amzn2023.0.5 and binutils
   2.41-50.amzn2023.0.5.
3. The resolver chose `kernel6.18-headers` for the kernel headers. With those headers `ms_check`
   matched its pin but the three ms2 binaries did not; in `ms2-bench` the only differing code was
   the `PerfCounters` constructor (four instructions and its unwind entry, plus the build ID), because
   `struct perf_event_attr` is 136 bytes in the 6.18 headers and 128 bytes in the 6.1 ones. The
   freeze box ran kernel 6.1.189-236.400, so a copy of the sandbox had the headers replaced with
   `kernel-headers-6.1.189-236.400.amzn2023`, and the build ran again in that copy (`run_ms2.sh`).
4. `run_ms2.sh` bound the unpacked freeze code tarball (`results/code/freeze-code.before-fig3count-fix.tar.gz`,
   which ran the build job) at `/work/code` and the three source bundles at their staged paths, and
   ran the freeze job's command line: `build_tools.py --parts ms2 --ms2-repo ... --ms2-commit
   2885ee86... --root /work/stage/work/freeze-2885ee8-3c8c8967/build-tools/outputs --work ...`.

Peak memory was 0.9 GB.
