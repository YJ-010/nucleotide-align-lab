"""Add public spliced-transcript annotations without duplicating verified CDS.

Usage: python build_transcripts.py config.json CDS_DATA_DIR OUTPUT_DIR
Only the exact reference transcript IDs already in the CDS catalog are included.
The existing CDS files and source FASTA/GFF files are read-only.
"""
import collections
import gzip
import hashlib
import json
import pathlib
import sys
from build_reference import attrs, digest, fasta, op, rc, shard, writegz


def unavailable(reason):
    return {'status': 'unavailable', 'reason': reason}


def assemble(record, genome, exons, utrs, transcript):
    """Return prefix/suffix only after coordinate-based CDS validation."""
    cds_parts = record['cdsSegments']
    chrom, strand = record['chrom'], record['strand']
    if exons:
        segments = sorted(set((p[1], p[2]) for p in exons))
        source = 'exons'
        if any(a[1] >= b[0] for a, b in zip(segments, segments[1:])):
            return unavailable('overlapping_exons')
        features = exons
    elif utrs:
        if any(max(s, u[1]) <= min(e, u[2]) for s, e, _ in cds_parts for u in utrs):
            return unavailable('utr_cds_overlap')
        segments = sorted(set([(s, e) for s, e, _ in cds_parts] + [(p[1], p[2]) for p in utrs]))
        source = 'cds_and_utr'
        features = utrs
    else:
        return unavailable('no_exon_or_utr')
    if chrom not in genome or strand not in ('+', '-'):
        return unavailable('missing_chromosome_or_strand')
    if any(p[0] != chrom or p[3] != strand for p in features):
        return unavailable('inconsistent_chromosome_or_strand')
    if any(s < 1 or e > len(genome[chrom]) or s > e for s, e in segments):
        return unavailable('invalid_coordinates')
    if any(s < transcript['start'] or e > transcript['end'] for s, e in segments):
        return unavailable('outside_transcript')
    # Union adjacent CDS/UTR intervals; gaps are introns, never filled with genomic DNA.
    merged = []
    for s, e in segments:
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    oriented = merged if strand == '+' else list(reversed(merged))
    sequence = ''.join(genome[chrom][s-1:e] if strand == '+' else rc(genome[chrom][s-1:e]) for s, e in oriented)
    if not sequence or any(c not in 'ACGTRYSWKMBDHVN' for c in sequence):
        return unavailable('invalid_sequence')
    coding = cds_parts if strand == '+' else list(reversed(cds_parts))
    positions = []
    for cs, ce, phase in coding:
        cursor = 0
        found = False
        for s, e in oriented:
            if s <= cs <= ce <= e:
                start = cursor + (cs-s if strand == '+' else e-ce)
                positions.append([start, start+ce-cs])
                found = True
                break
            cursor += e-s+1
        if not found:
            return unavailable('cds_outside_exon')
    if any(a[1]+1 != b[0] for a, b in zip(positions, positions[1:])):
        return unavailable('noncontiguous_cds_in_transcript')
    start, end = positions[0][0]+coding[0][2], positions[-1][1]+1
    if sequence[start:end] != record['cds']:
        return unavailable('cds_sequence_mismatch')
    prefix, suffix = sequence[:start], sequence[end:]
    warnings = []
    if not prefix: warnings.append('no_5utr')
    if not suffix: warnings.append('no_3utr')
    if coding[0][2]: warnings.append('partial_cds_phase')
    if any(c not in 'ACGT' for c in sequence): warnings.append('ambiguous_bases')
    return {'status': 'available', 'method': source, 'prefix': prefix, 'suffix': suffix,
            'length': len(sequence), 'cdsLength': len(record['cds']), 'cdsStart': start+1,
            'cdsEnd': end, 'exons': merged, 'warnings': warnings}


def build(config, base, out):
    base, out = pathlib.Path(base), pathlib.Path(out)
    out.mkdir(parents=True, exist_ok=True)
    original = json.loads((base/'manifest.json').read_text())
    known = {v['id']: v for v in original['versions']}
    report = {'schema': 1, 'built': '2026-10-19', 'alphabet': 'DNA (T), spliced mRNA/cDNA representation',
              'storage': 'verified CDS plus annotated transcript prefix and suffix', 'versions': []}
    publish = []
    for conf in config:
        version = conf['id']
        for key in ['genome', 'gff']:
            if digest(conf[key]) != known[version]['sha256'][key]:
                raise ValueError('Reference source changed: '+version+' '+key)
        print('TRANSCRIPTS', version, flush=True)
        tx, exons, utrs = {}, collections.defaultdict(list), collections.defaultdict(list)
        with op(conf['gff']) as fh:
            for line in fh:
                if line.startswith('#'): continue
                a = line.rstrip().split('\t')
                if len(a) != 9: continue
                ch, _, kind, s, e, _, st, _, at = a
                at = attrs(at)
                if kind in ('mRNA', 'transcript'):
                    tx[at['ID']] = {'start': int(s), 'end': int(e)}
                elif kind in ('exon', 'five_prime_UTR', 'three_prime_UTR', 'UTR'):
                    for parent in at.get('Parent', '').split(','):
                        (exons if kind == 'exon' else utrs)[parent].append((ch, int(s), int(e), st))
        genome = fasta(conf['genome']) if exons or utrs else {}
        counts, reasons = collections.Counter(), collections.Counter()
        all_rows = {}
        for path in sorted((base/version).glob('*.json.gz')):
            if path.name == 'index.json.gz': continue
            records = json.loads(gzip.decompress(path.read_bytes()))
            for tid, record in records.items():
                result = assemble(record, genome, exons.get(tid, []), utrs.get(tid, []), tx.get(tid, {}))
                all_rows[tid] = result
                counts[result['status']] += 1
                if result['status'] == 'available':
                    counts[result['method']] += 1
                    counts['with_5prime_extension'] += bool(result['prefix'])
                    counts['with_3prime_extension'] += bool(result['suffix'])
                    counts['no_annotated_utr'] += not result['prefix'] and not result['suffix']
                    # Storage must reconstruct the entire independently extracted transcript exactly.
                    seq = result['prefix']+record['cds']+result['suffix']
                    assert len(seq) == result['length']
                    assert seq[result['cdsStart']-1:result['cdsEnd']] == record['cds']
                else:
                    reasons[result['reason']] += 1
        assert len(all_rows) == known[version]['transcriptCount']
        count = 8
        target = out/version
        target.mkdir(exist_ok=True)
        while True:
            buckets = [dict() for _ in range(count)]
            for tid, row in all_rows.items(): buckets[int(shard(tid), 16) % count][tid] = row
            encoded = [gzip.compress(json.dumps(b, separators=(',', ':')).encode(), 6, mtime=0) for b in buckets]
            if max(map(len, encoded)) <= 1500000 or count == 64: break
            count *= 2
        for i, content in enumerate(encoded):
            path = target/f'{i:02x}.json.gz'
            path.write_bytes(content)
            publish.append({'path': str(path.relative_to(out)), 'bytes': len(content),
                            'sha256': hashlib.sha256(content).hexdigest(),
                            'gitsha': hashlib.sha1(b'blob '+str(len(content)).encode()+b'\0'+content).hexdigest()})
        meta = {'id': version, 'shards': count, 'counts': dict(counts), 'unavailableReasons': dict(reasons),
                'source': conf['source'], 'sha256': known[version]['sha256']}
        report['versions'].append(meta)
        print(json.dumps(meta), flush=True)
    report['files'] = len(publish)
    report['bytes'] = sum(p['bytes'] for p in publish)
    (out/'manifest.json').write_text(json.dumps(report, indent=2)+'\n')
    (out/'publish-verified.json').write_text(json.dumps(publish, indent=2)+'\n')
    print('DONE', report['files'], report['bytes'], flush=True)
    return report


if __name__ == '__main__':
    build(json.load(open(sys.argv[1])), sys.argv[2], sys.argv[3])
