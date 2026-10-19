"""Build public, versioned gene/CDS/upstream shards and reciprocal protein hits.
Usage: python build_reference.py config.json OUTPUT_DIR [--threads 4]
Config: [{id,label,genome,gff,source,cds(optional),protein(optional)}]. No private data needed.
Requires DIAMOND on PATH. Python standard library only. Raw references are read-only.
"""
import argparse, concurrent.futures, collections, gzip, hashlib, json, os, pathlib, subprocess, time, urllib.parse
COMP=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
def rc(s):return s.translate(COMP)[::-1]
def op(p):return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)
def fasta(p):
 out={};name=None;seq=[]
 with op(p) as f:
  for l in f:
   if l.startswith('>'):
    if name:out[name]=''.join(seq).upper()
    name=l[1:].split()[0];seq=[]
   else:seq.append(l.strip())
 if name:out[name]=''.join(seq).upper()
 return out
def attrs(s):
 return {k.strip():urllib.parse.unquote(v.strip()) for x in s.split(';') if '=' in x for k,v in [x.split('=',1)]}
CODONS=dict(zip([a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'], 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
def protein(s):return ''.join(CODONS.get(s[i:i+3],'X') for i in range(0,len(s)-2,3)).rstrip('*')
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def shard(key):
 h=2166136261
 for ch in key.lower():h=((h^ord(ch))*16777619)&0xffffffff
 return f'{h%64:02x}'
def writegz(p,obj):
 b=json.dumps(obj,separators=(',',':'),ensure_ascii=False).encode();p.write_bytes(gzip.compress(b,6,mtime=0))
def extract(c,out):
 print('EXTRACT',c['id'],flush=True)
 genome=fasta(c['genome']);genes={};tx={};parts=collections.defaultdict(list);protein_keys={}
 with op(c['gff']) as f:
  for l in f:
   if l.startswith('#'):continue
   a=l.rstrip().split('\t')
   if len(a)!=9:continue
   chrom,_,kind,s,e,_,strand,phase,at=a;at=attrs(at)
   if kind=='gene':genes[at['ID']]={'chrom':chrom,'start':int(s),'end':int(e),'strand':strand,'attrs':at}
   elif kind in ('mRNA','transcript'):tx[at['ID']]={'gene':at.get('Parent','').split(',')[0],'chrom':chrom,'start':int(s),'end':int(e),'strand':strand,'attrs':at}
   elif kind=='CDS':
    for parent in at.get('Parent','').split(','):
     parts[parent].append((chrom,int(s),int(e),strand,int(phase) if phase.isdigit() else 0));protein_keys[parent]=at.get('protein_id',parent)
 published_proteins=fasta(c['protein']) if c.get('protein') else {}
 supplied=fasta(c['cds']) if c.get('cds') else {}; records={};bygene=collections.defaultdict(list);skips=collections.Counter();audit=collections.Counter()
 for tid,t in tx.items():
  ps=sorted(parts.get(tid,[]),key=lambda x:x[1]);g=genes.get(t['gene']);ch=t['chrom'];st=t['strand']
  if not g or not ps:skips['no_gene_or_cds']+=1;continue
  if ch not in genome or st not in ('+','-'):skips['missing_chromosome_or_strand']+=1;continue
  if any(p[0]!=ch or p[3]!=st or p[1]<1 or p[2]>len(genome[ch]) or p[1]>p[2] for p in ps):skips['invalid_coordinates']+=1;continue
  if any(a[2]>=b[1] for a,b in zip(ps,ps[1:])):skips['overlapping_CDS']+=1;continue
  seq=''.join(genome[ch][p[1]-1:p[2]] for p in ps)
  if st=='-':seq=rc(seq)
  first=ps[0] if st=='+' else ps[-1];phase=first[4];seq=seq[phase:]
  canonical=g['attrs'].get('locus_tag') or g['attrs'].get('Name') or t['gene']
  if c['id']=='red5':canonical=canonical.removeprefix('CEY00_')
  if not canonical.startswith(('Ach','Actinidia','DTZ79_','Acc')):canonical=t['gene'].removeprefix('gene-')
  aa=protein(seq);warnings=[]
  if phase:warnings.append('CDS 首段 phase 非零，可能是不完整模型')
  if not seq.startswith('ATG'):warnings.append('CDS 起点不是 ATG；上游区以注释 CDS 边界为基准')
  if len(seq)%3:warnings.append('CDS 长度不是 3 的倍数')
  if seq[-3:] not in ('TAA','TAG','TGA'):warnings.append('CDS 无标准终止密码子')
  if '*' in aa:warnings.append('CDS 翻译含内部终止密码子')
  if any(x not in 'ACGT' for x in seq):warnings.append('CDS 含不确定碱基')
  if st=='+':anchor=ps[0][1]-1;pstart=max(0,anchor-2000);pend=anchor;up=genome[ch][pstart:pend]
  else:anchor=ps[-1][2];pstart=anchor;pend=min(len(genome[ch]),anchor+2000);up=rc(genome[ch][pstart:pend])
  if len(up)<2000:warnings.append('靠近序列边界，上游区不足 2000 bp')
  if any(p[1]<t['start'] or p[2]>t['end'] for p in ps):warnings.append('源注释 CDS 超出转录本边界')
  match=supplied.get(tid) or supplied.get(canonical)
  if match is not None:
   audit['supplied_match' if match==seq else 'supplied_mismatch']+=1
   if match!=seq:warnings.append('基因组提取 CDS 与配套 CDS 文件不完全一致')
  published=published_proteins.get(protein_keys.get(tid,tid)) or published_proteins.get(canonical)
  if published is not None:
   agrees=aa==published.rstrip('*');audit['publishedProteinMatch' if agrees else 'publishedProteinMismatch']+=1
   if not agrees:warnings.append('基因组提取 CDS 的翻译与发布蛋白不完全一致，请核对注释模型')
  rec={'id':canonical,'transcript':tid,'chrom':ch,'strand':st,'geneStart':g['start'],'geneEnd':g['end'],'cdsSegments':[[p[1],p[2],p[4]] for p in ps], 'upstreamStart':pstart+1,'upstreamEnd':pend,'cds':seq,'upstream':up,'warnings':warnings,'transcriptLabel':t['attrs'].get('orig_protein_id',tid).split('|')[-1], 'description':t['attrs'].get('product',''), 'aliases':list(set([canonical,g['attrs'].get('locus_tag',''),t['gene'],tid,t['attrs'].get('Name',''),t['attrs'].get('orig_protein_id','').split('|')[-1],t['attrs'].get('orig_transcript_id','').split('|')[-1]]))}
  bygene[canonical].append((len(seq),tid,rec,aa))
 # Default gene query uses longest CDS; transcript ID queries always preserve its actual sequence.
 index={};pep=[];representatives={}
 for gid,items in bygene.items():
  items.sort(key=lambda x:(-x[0],x[1]));representatives[gid]=items[0][1]
  for _,tid,r,aa in items:
   r['isoforms']=[x[1] for x in items];r['representative']=items[0][1];records[tid]=r
   for alias in r['aliases']:
    if alias and alias != gid:index[alias.lower()]=tid
   index[tid.lower()]=tid
  best=items[0];index[gid.lower()]=best[1]
  if best[3] and len(best[3])>=20 and '*' not in best[3]:pep.append('>'+gid+'\n'+best[3]+'\n')
 # Resolve raw gene IDs to the selected representative, not an arbitrary isoform.
 for tid,r in records.items():
  for alias in r['aliases']:
   if alias in genes:index[alias.lower()]=representatives[r['id']]
 target=out/c['id'];target.mkdir(parents=True,exist_ok=True)
 (target/'proteins.fa').write_text(''.join(pep));(target/'records.json').write_text(json.dumps(records,separators=(',',':')))
 writegz(target/'index.json.gz',index)
 meta={k:c[k] for k in ['id','label','source']};meta.update({'geneCount':len(bygene),'transcriptCount':len(records),'proteinCount':len(pep),'skippedModels':dict(skips),'sequenceAudit':dict(audit),'sha256':{'genome':digest(c['genome']),'gff':digest(c['gff'])}})
 (target/'meta.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2))
 print('EXTRACTED',meta,flush=True);return meta

def hits(task,out,threads):
 q,t=task;dest=out/f'{q}__{t}.tsv'
 if not dest.exists():
  cmd=['diamond','blastp','--query',str(out/q/'proteins.fa'),'--db',str(out/t/'protein_db'),'--out',str(dest),'--threads',str(threads),'--sensitive','--evalue','1e-10','--max-target-seqs','5','--max-hsps','1','--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore','--quiet']
  subprocess.run(cmd,check=True)
 print('HITS',q,t,dest.stat().st_size,flush=True)
 return q,t

def main():
 ap=argparse.ArgumentParser();ap.add_argument('config');ap.add_argument('out');ap.add_argument('--threads',type=int,default=4);args=ap.parse_args();conf=json.load(open(args.config));out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True);metas=[]
 for c in conf:
  p=out/c['id']/'meta.json'
  meta=json.loads(p.read_text()) if p.exists() else extract(c,out);metas.append(meta)
  db=out/c['id']/'protein_db'
  if not pathlib.Path(str(db)+'.dmnd').exists():subprocess.run(['diamond','makedb','--in',str(out/c['id']/'proteins.fa'),'--db',str(db),'--quiet'],check=True)
 versions=[c['id'] for c in conf]
 with concurrent.futures.ThreadPoolExecutor(3) as pool:list(pool.map(lambda x:hits(x,out,args.threads),[(q,t) for q in versions for t in versions if q!=t]))
 allhits={};best={}
 for q in versions:
  for t in versions:
   if q==t:continue
   h=collections.defaultdict(list)
   for line in (out/f'{q}__{t}.tsv').read_text().splitlines():
    a=line.split('\t');qc=(int(a[7])-int(a[6])+1)/int(a[4])*100;tc=(int(a[9])-int(a[8])+1)/int(a[5])*100
    if qc<40 or tc<40 or float(a[2])<30:continue
    h[a[0]].append({'id':a[1],'identity':float(a[2]),'queryCoverage':round(qc,2),'targetCoverage':round(tc,2),'evalue':float(a[10]),'bitscore':float(a[11])})
   for gid,arr in h.items():
    arr.sort(key=lambda x:(-x['bitscore'],x['evalue'],x['id']));best[(q,t,gid)]={x['id'] for x in arr if abs(x['bitscore']-arr[0]['bitscore'])<1e-6}
   allhits[(q,t)]=h
 for q in versions:
  records=json.loads((out/q/'records.json').read_text());buckets=collections.defaultdict(dict)
  for tid,r in records.items():
   r['homologs']={}
   for t in versions:
    if t==q:continue
    arr=allhits[(q,t)].get(r['id'],[])[:3]
    for h in arr:
     forward=best.get((q,t,r['id']),set());back=best.get((t,q,h['id']),set())
     h['reciprocal']=r['id'] in back and h['id'] in forward
     h['tiedBest']=h['reciprocal'] and (len(forward)>1 or len(back)>1)
    r['homologs'][t]=arr
   buckets[shard(tid)][tid]=r
  for key,rows in buckets.items():writegz(out/q/(key+'.json.gz'),rows)
 manifest={'schema':1,'built':time.strftime('%Y-%m-%d'),'upstreamLimit':2000,'upstreamAnchor':'annotated CDS 5-prime boundary (ATG where complete), excludes CDS; transcriptional orientation','homology':'DIAMOND blastp sensitive; E <= 1e-10; identity >= 30%; query/target coverage >= 40%; top 3 per reference; reciprocal top hits are candidates, not proven orthologs','representative':'longest genomic CDS per gene; exact transcript identifiers retain isoform; homolog search uses representative proteins','versions':metas}
 (out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
 export=[p for p in out.rglob('*') if p.is_file() and (p.name.endswith('.json.gz') or p.name=='manifest.json')]
 (out/'publish-list.json').write_text(json.dumps([{'path':str(p.relative_to(out)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in export]))
 print('DONE',len(export),sum(p.stat().st_size for p in export),flush=True)
if __name__=='__main__':main()
