"""Add verified, public genomic flanks around existing promoter targets.

Usage: python build_primer_flanks.py CONFIG ORIGINAL_DATA EXPECTED_SHAS OUTPUT
Original FASTA, annotations and shards are read-only. Output must be new.
"""
import argparse, base64, collections, gzip, hashlib, json, pathlib

def revcomp(seq):
    return seq.translate(str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN'))[::-1]

def shard(key):
    value=2166136261
    for char in key.lower():
        value=((value^ord(char))*16777619)&0xffffffff
    return value%16

def genome_file(path):
    opener=gzip.open if path.endswith('.gz') else open
    seqs={}; name=None; parts=[]
    with opener(path,'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None: seqs[name]=''.join(parts).upper()
                name=line[1:].split()[0]; parts=[]
            else: parts.append(line.strip())
    if name is not None: seqs[name]=''.join(parts).upper()
    return seqs

def checksum(path):
    value=hashlib.sha256()
    with open(path,'rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):value.update(block)
    return value.hexdigest()

def run(config_path,data_path,expected_path,output_path):
    configs=json.loads(pathlib.Path(config_path).read_text())
    data=pathlib.Path(data_path); out=pathlib.Path(output_path)
    out.mkdir(parents=True,exist_ok=False)
    expected=json.loads(pathlib.Path(expected_path).read_text())
    manifest=json.loads((data/'manifest.json').read_text())
    hashes={r['id']:r['sha256']['genome'] for r in manifest['versions']}
    report={'schema':1,'built':'2026-10-26','flankLimit':200,'shards':16,
            'coordinateSystem':'1-based inclusive; all sequences in gene transcriptional orientation',
            'sources':[],'files':[]}
    samples=[]
    for config in configs:
        version=config['id']; path=config['genome']
        observed=checksum(path)
        assert observed==hashes[version],f'{version}: genome checksum mismatch'
        genome=genome_file(path); groups=[{} for _ in range(16)]
        counts=collections.Counter()
        for source_file in sorted((data/version).glob('??.json.gz')):
            raw=source_file.read_bytes()
            gitsha=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
            relative=f'{version}/{source_file.name}'
            assert expected[relative]==gitsha,f'{relative}: differs from published source'
            records=json.loads(gzip.decompress(raw))
            for transcript,record in records.items():
                chrom=genome[record['chrom']]; start=record['upstreamStart']-1; end=record['upstreamEnd']
                target=chrom[start:end]
                low=chrom[max(0,start-200):start]; high=chrom[end:min(len(chrom),end+200)]
                if record['strand']=='-':target=revcomp(target);left,right=revcomp(high),revcomp(low)
                else:left,right=low,high
                assert target==record['upstream'],f'{version}/{transcript}: promoter mismatch'
                entry={'l':left,'r':right,'s':start+1,'e':end,'n':len(target),'c':record['chrom'],'st':record['strand']}
                groups[shard(transcript)][transcript]=entry
                counts['transcripts']+=1;counts['full_200bp_flanks']+=int(len(left)==200 and len(right)==200)
                counts['sufficient_18bp_flanks']+=int(len(left)>=18 and len(right)>=18)
                counts['negative_strand']+=int(record['strand']=='-')
                if record['id']=='Actinidia03760' or (len(samples)<6 and len(left)==200 and record['strand']=='-'):
                    samples.append({'version':version,'record':record,'promoterFlanks':entry})
        (out/version).mkdir()
        for i,rows in enumerate(groups):
            encoded=json.dumps(rows,separators=(',',':'),ensure_ascii=False).encode()
            payload=gzip.compress(encoded,compresslevel=9,mtime=0)
            relative=f'{version}/{i:02x}.json.gz';(out/relative).write_bytes(payload)
            (out/(relative+'.b64')).write_text(base64.b64encode(payload).decode())
            report['files'].append({'path':relative,'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest(),
                                    'gitsha':hashlib.sha1(b'blob '+str(len(payload)).encode()+b'\0'+payload).hexdigest()})
        report['sources'].append({'id':version,'source':config['source'],'genomeSha256':observed,**counts})
        print(json.dumps({'version':version,**counts}),flush=True)
    (out/'manifest.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    (out/'qa-samples.json').write_text(json.dumps(samples,ensure_ascii=False))
    print(json.dumps({'status':'passed','files':len(report['files']),'bytes':sum(r['bytes'] for r in report['files'])}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');p.add_argument('data');p.add_argument('expected');p.add_argument('output')
    a=p.parse_args();run(a.config,a.data,a.expected,a.output)
