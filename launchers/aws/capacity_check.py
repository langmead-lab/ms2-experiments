#!/usr/bin/env python3
"""Read-only check of spot price, spot capacity and vCPU quota for the freeze fleet (D4).

Usage:
  capacity_check.py [--config config.json] [--boxes N] [--days 7] [--profile PROFILE]

Prints, for each instance type of config.json in its region:
- the spot price per availability zone over the last --days (minimum, median, latest);
- the spot placement score (1 to 10) for --boxes boxes of that type, per zone and for the region;
- the account's spot vCPU quota for standard instances (L-34B43A08) and what --boxes needs.

Nothing is launched or changed. Needs credentials that allow ec2:DescribeSpotPriceHistory,
ec2:GetSpotPlacementScores and servicequotas:GetServiceQuota.
"""
import argparse
from collections import defaultdict
import datetime as dt
import json
from pathlib import Path
import statistics
import subprocess
import sys

HERE = Path(__file__).resolve().parent
SPOT_QUOTA = 'L-34B43A08'  # All Standard (A, C, D, H, I, M, R, T, Z) Spot Instance Requests, in vCPUs


def aws(profile, region, *args):
    argv = ['aws', *(['--profile', profile] if profile else []), '--region', region, '--output', 'json', *args]
    result = subprocess.run(argv, capture_output=True, text=True)
    if result.returncode:
        raise SystemExit(f"{' '.join(args[:2])} failed: {result.stderr.strip()}")
    return json.loads(result.stdout or '{}')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', default=str(HERE/'config.json'))
    p.add_argument('--boxes', type=int, help='target box count (default config box_count)')
    p.add_argument('--days', type=int, default=7)
    p.add_argument('--profile')
    args = p.parse_args()
    config = json.loads(Path(args.config).read_text())
    region, types = config['region'], config['instance_types']
    boxes = args.boxes or config['box_count']
    start = (dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=args.days)).strftime('%Y-%m-%dT%H:%M:%SZ')
    vcpus = {}
    for itype in types:
        info = aws(args.profile, region, 'ec2', 'describe-instance-types', '--instance-types', itype)
        vcpus[itype] = info['InstanceTypes'][0]['VCpuInfo']['DefaultVCpus']
        history = aws(args.profile, region, 'ec2', 'describe-spot-price-history', '--instance-types', itype,
                      '--product-descriptions', 'Linux/UNIX', '--start-time', start)['SpotPriceHistory']
        by_zone = defaultdict(list)
        for row in sorted(history, key=lambda r: r['Timestamp']):
            by_zone[row['AvailabilityZone']].append(float(row['SpotPrice']))
        print(f'{itype} ({vcpus[itype]} vCPUs), spot $/h over {args.days} days')
        for zone, prices in sorted(by_zone.items()):
            print(f'  {zone}: min {min(prices):.3f}  median {statistics.median(prices):.3f}  latest {prices[-1]:.3f}')
        if not by_zone:
            print('  no spot price history')
        for single in (False, True):
            scores = aws(args.profile, region, 'ec2', 'get-spot-placement-scores', '--instance-types', itype,
                         '--target-capacity', str(boxes), '--target-capacity-unit-type', 'units',
                         '--region-names', region, *(['--single-availability-zone'] if single else []))
            for row in scores.get('SpotPlacementScores', []):
                where = row.get('AvailabilityZoneId', row['Region'])
                print(f"  placement score for {boxes} boxes in {where}: {row['Score']}")
    quota = aws(args.profile, region, 'service-quotas', 'get-service-quota', '--service-code', 'ec2',
                '--quota-code', SPOT_QUOTA)['Quota']['Value']
    need = boxes * max(vcpus.values())
    print(f'spot vCPU quota (standard): {quota:.0f}; {boxes} boxes need {need}'
          + ('' if quota >= need else f'; at most {int(quota // max(vcpus.values()))} boxes fit'))
    return 0 if quota >= need else 1


if __name__ == '__main__':
    sys.exit(main())
