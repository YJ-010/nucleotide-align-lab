import { reverseComplement } from './align-core.mjs?v=20261019-clone1';
export { reverseComplement };
export const VERSIONS=[['hy2','Hong Yang v2'],['hy3','Hong Yang v3'],['red5','Red5'],['white','White v1.0']];
export function shardFor(id){let h=2166136261;for(const c of id.toLowerCase())h=Math.imul(h^c.charCodeAt(0),16777619)>>>0;return (h%64).toString(16).padStart(2,'0');}
export function promoterSequence(record,length=2000){
 if(!Number.isInteger(length)||length<60||length>2000)throw new Error('上游克隆片段长度请填 60–2000 bp。');
 if(!record.upstream)throw new Error('该记录未提供上游区序列。');
 const sequence=record.upstream.slice(-length), n=sequence.length;
 const located=Number.isInteger(record.upstreamStart)&&Number.isInteger(record.upstreamEnd);
 const start=located?(record.strand==='+'?record.upstreamEnd-n+1:record.upstreamStart):null;
 const end=located?(record.strand==='+'?record.upstreamEnd:record.upstreamStart+n-1):null;
 return {sequence,start,end,requested:length,truncated:n<length};
}
export function cloningTemplate(record,kind,length=2000,removeStop=false){
 if(kind==='promoter')return {...promoterSequence(record,length),kind};
 if(!record.cds)throw new Error('该记录未提供 CDS 序列。');
 let sequence=record.cds,removedStop=false;
 if(removeStop){if(!/T(?:AA|AG|GA)$/.test(sequence))throw new Error('该 CDS 没有标准终止密码子，无法自动去除。');sequence=sequence.slice(0,-3);removedStop=true;}
 return {sequence,kind:'cds',removedStop};
}
// End-to-end Needleman–Wunsch. Exact canonical matches only; gaps and ambiguities stay in denominator.
export function compareDNA(a,b){
 if(!/^[ACGTRYSWKMBDHVN]+$/.test(a)||!/^[ACGTRYSWKMBDHVN]+$/.test(b))throw new Error('比较序列为空或含非法字符。');
 if(a.length*b.length>36000000||Math.max(a.length,b.length)>12000)throw new Error('该 CDS 超过单次网页比较上限，请下载 FASTA 后使用专业软件比较。');
 const n=a.length,m=b.length,w=m+1,trace=new Uint8Array((n+1)*w);let prev=new Int32Array(w),cur=new Int32Array(w);
 for(let j=1;j<=m;j++){prev[j]=-2*j;trace[j]=2;}
 for(let i=1;i<=n;i++){
  cur[0]=-2*i;trace[i*w]=1;
  for(let j=1;j<=m;j++){
   const exact=a[i-1]===b[j-1]&&'ACGT'.includes(a[i-1]);
   const d=prev[j-1]+(exact?2:-1),u=prev[j]-2,l=cur[j-1]-2;
   if(d>=u&&d>=l){cur[j]=d;}else if(u>=l){cur[j]=u;trace[i*w+j]=1;}else{cur[j]=l;trace[i*w+j]=2;}
  }[prev,cur]=[cur,prev];
 }
 let i=n,j=m,ar=[],br=[],matches=0,paired=0,gaps=0,mismatches=0,ambiguous=0;
 while(i||j){const t=trace[i*w+j];if(i&&j&&t===0){const x=a[--i],y=b[--j];ar.push(x);br.push(y);paired++;if(x===y&&'ACGT'.includes(x))matches++;else if(!'ACGT'.includes(x)||!'ACGT'.includes(y))ambiguous++;else mismatches++;}else if(i&&(j===0||t===1)){ar.push(a[--i]);br.push('-');gaps++;}else{ar.push('-');br.push(b[--j]);gaps++;}}
 const columns=ar.length;return {alignedA:ar.reverse().join(''),alignedB:br.reverse().join(''),identity:100*matches/columns,pairedIdentity:paired?100*matches/paired:0,coverageA:100*paired/n,coverageB:100*paired/m,matches,paired,gaps,mismatches,ambiguous,columns};
}
// SantaLucia unified nearest-neighbor parameters; salt correction to entropy.
// Sodium-equivalent correction for Mg/dNTP: von Ahsen 2001. Cross-checked against Biopython Tm_NN, DNA_NN3, saltcorr=5.
const NN={AA:[-7.9,-22.2],TT:[-7.9,-22.2],AT:[-7.2,-20.4],TA:[-7.2,-21.3],CA:[-8.5,-22.7],TG:[-8.5,-22.7],GT:[-8.4,-22.4],AC:[-8.4,-22.4],CT:[-7.8,-21],AG:[-7.8,-21],GA:[-8.2,-22.2],TC:[-8.2,-22.2],CG:[-10.6,-27.2],GC:[-9.8,-24.4],GG:[-8,-19.9],CC:[-8,-19.9]};
export function meltingTemp(seq,{salt=50,mg=1.5,dntp=.2,concentration=250}={}){
 if(!/^[ACGT]{2,80}$/.test(seq))throw new Error('引物 Tm 仅支持 2–80 个明确 DNA 碱基。');
 if(![salt,mg,dntp,concentration].every(Number.isFinite)||salt<=0||mg<0||dntp<0||concentration<=0)throw new Error('盐浓度和引物浓度参数无效。');
 let h=0,s=0;
 for(let i=0;i<seq.length-1;i++){const t=NN[seq.slice(i,i+2)];h+=t[0];s+=t[1];}
 for(const base of [seq[0],seq.at(-1)]){if('AT'.includes(base)){h+=2.3;s+=4.1;}else{h+=.1;s-=2.8;}}
 const symmetric=seq===reverseComplement(seq);if(symmetric)s-=1.4;
 const sodium=(salt+120*Math.sqrt(Math.max(0,mg-dntp)))/1000;
 s+=.368*(seq.length-1)*Math.log(sodium);
 const effective=concentration*1e-9/(symmetric?2:4);
 return 1000*h/(s+1.987*Math.log(effective))-273.15;
}
export function complementarity(a,b){
 const rc=reverseComplement(b);let any=0,end=0;
 for(let offset=-rc.length+1;offset<a.length;offset++){
  let run=0;for(let i=0;i<a.length;i++){
   const j=i-offset;if(j>=0&&j<rc.length&&a[i]===rc[j]){run++;any=Math.max(any,run);if(i===a.length-1||j-run+1===0)end=Math.max(end,run);}else run=0;
  }
 }
 return {any,end};
}
export function hairpinStem(s){let best=0;for(let i=0;i<s.length;i++)for(let j=i+6;j<s.length;j++){let k=0;while(i+k<j-k-3&&s[i+k]===reverseComplement(s[j-k]))k++;best=Math.max(best,k);}return best;}
function primerMetrics(seq,opts){const gc=100*[...seq].filter(c=>'GC'.includes(c)).length/seq.length;return {sequence:seq,length:seq.length,tm:meltingTemp(seq,opts),gc,homopolymer:Math.max(...(seq.match(/(.)\1*/g)||[]).map(x=>x.length)),hairpin:hairpinStem(seq),self:complementarity(seq,seq)};}
export function designPrimers(sequence,options={}){
 const target=options.target??60,minLen=options.minLen??18,maxLen=options.maxLen??30;
 if(!Number.isFinite(target)||target<50||target>70)throw new Error('目标 Tm 请设为 50–70 °C。');
 if(!Number.isInteger(minLen)||!Number.isInteger(maxLen)||minLen<18||maxLen>35||minLen>maxLen)throw new Error('引物长度范围须为 18–35 nt。');
 if(!/^[ACGTRYSWKMBDHVN]+$/.test(sequence)||sequence.length<60||sequence.length>12000)throw new Error('克隆模板须为 60–12,000 bp，且仅含 DNA 碱基。');
 const tails=[options.forwardTail||'',options.reverseTail||''].map(s=>s.replace(/\s/g,'').toUpperCase());
 if(tails.some(x=>!/^[ACGT]{0,60}$/.test(x)))throw new Error('5′ 附加序列只能含 A/C/G/T，每端最多 60 nt。');
 const lists=[[],[]];
 for(let side=0;side<2;side++)for(let len=minLen;len<=maxLen;len++){
  const seq=side===0?sequence.slice(0,len):reverseComplement(sequence.slice(-len));if(!/^[ACGT]+$/.test(seq))continue;
  const p=primerMetrics(seq,options),full=tails[side]+seq;
  p.fullSequence=full;p.tail=tails[side];p.hairpin=hairpinStem(full);p.self=complementarity(full,full);
  if(Math.abs(p.tm-target)>5||p.gc<25||p.gc>75||p.homopolymer>5||p.hairpin>=6||p.self.end>=5)continue;
  lists[side].push(p);
 }
 const pairs=[];
 for(const f of lists[0])for(const r of lists[1]){
  const delta=Math.abs(f.tm-r.tm),cross=complementarity(f.fullSequence,r.fullSequence);
  if(delta>3||cross.end>=5||f.length+r.length>sequence.length)continue;
  const score=Math.abs(f.tm-target)+Math.abs(r.tm-target)+delta+Math.abs(f.length-22)*.15+Math.abs(r.length-22)*.15+f.self.end+r.self.end+cross.end;
  pairs.push({forward:f,reverse:r,tmDifference:delta,heterodimer:cross,productSize:sequence.length,productWithTails:sequence.length+tails[0].length+tails[1].length,score});
 }
 pairs.sort((a,b)=>a.score-b.score||a.forward.length-b.forward.length);return {pairs:pairs.slice(0,3),parameters:{target,minLen,maxLen,salt:options.salt??50,mg:options.mg??1.5,dntp:options.dntp??.2,concentration:options.concentration??250},templateAmbiguous:(sequence.match(/[^ACGT]/g)||[]).length,reason:pairs.length?'':'固定片段边界下未找到同时满足 Tm、GC、长度及连续互补初筛的引物。请调整目标 Tm／长度，或更改克隆片段边界；不会自动截短您的目标片段。'};
}
export function bindingCheck(primer,sequence,reverse=false){
 const target=reverse?reverseComplement(sequence):sequence;let best=null;
 for(let i=0;i<=target.length-primer.length;i++){
  let mismatches=0,endMismatches=0;
  for(let j=0;j<primer.length;j++){if(primer[j]!==target[i+j]||!'ACGT'.includes(target[i+j])){mismatches++;if(j>=primer.length-5)endMismatches++;}}
  const score=endMismatches*10+mismatches;
  if(!best||score<best.score){best={score,mismatches,endMismatches,matches:1,start:reverse?sequence.length-i-primer.length+1:i+1,end:reverse?sequence.length-i:i+primer.length};}
  else if(score===best.score)best.matches++;
 }return best;
}
export function fastaText(rows,kind,length){return rows.map(r=>{const t=cloningTemplate(r,kind,length);return `>${r.version}|${r.id}|${r.transcript}|${kind}${kind==='promoter'&&t.start!==null?`|${r.chrom}:${t.start}-${t.end}(${r.strand})`:r.manual?'|user_provided':''}\n${t.sequence.match(/.{1,80}/g)?.join('\n')||''}`;}).join('\n');}
