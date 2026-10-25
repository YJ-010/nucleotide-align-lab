import {cloningTemplate,designPrimers,meltingTemp,complementarity,hairpinStem,reverseComplement,hasTandemRepeat} from './clone-core.mjs?v=20261026-boundary2';

// Coordinates below are zero-based half-open on a 5′→3′ cDNA or genomic context.
// A primer's reported start/end always describe its annealing site on that context.
export function cloningContext(record,kind='cds',length=2000,override=''){
 const template=cloningTemplate(record,kind,length,false),target=template.sequence;
 if(kind==='cds'&&(!target.startsWith('ATG')||!/(TAA|TAG|TGA)$/.test(target)||target.length%3))throw new Error('该 CDS 的 ATG、末端终止密码子或三联体长度不完整，请核对转录本后再设计完整 CDS 引物。');
 if(kind==='promoter'&&template.truncated)throw new Error(`该记录受参考序列边界限制，只有 ${target.length} bp，无法覆盖所要求的完整 ${length} bp 启动子候选区。`);
 const base={target,kind,template,sequence:target,start:0,end:target.length,source:'',reason:''};
 if(override){
  if(!/^[ACGTRYSWKMBDHVN]+$/.test(override)||override.length>30000)throw new Error('补充序列须为不超过 30,000 bp 的一条 DNA 序列。');
  const start=override.indexOf(target);
  if(start<0||start!==override.lastIndexOf(target))throw new Error('补充序列必须按 5′→3′ 方向且仅一次完整包含当前目标，不能缺失、替换或反向输入目标。');
  return {...base,sequence:override,start,end:start+target.length,source:kind==='cds'?'用户提供的 cDNA 两侧序列（来源待核对）':'用户提供的连续基因组两侧序列（来源待核对）'};
 }
 if(kind==='cds'){
  const m=record.mrnaInfo,start=(m?.cdsStart??0)-1,end=m?.cdsEnd;
  if(!record.mrna||start<0||!Number.isInteger(end)||record.mrna.slice(start,end)!==target)return {...base,reason:'没有可核对的两侧 UTR。可补充含完整 CDS 及两侧 UTR 的 cDNA 序列；不会用含内含子的基因组侧翼代替。'};
  return {...base,sequence:record.mrna,start,end,source:record.manual?'用户提供的 mRNA / cDNA（来源待核对）':'同一转录本的已注释 5′UTR 与 3′UTR（cDNA）'};
 }
 if(kind==='promoter'){
  const f=record.promoterFlanks;
  if(!f)return {...base,reason:record.flankError||'暂无启动子外侧数据。可补充含完整目标及两侧序列的连续基因组片段。'};
  if(f.n!==record.upstream.length||f.s!==record.upstreamStart||f.e!==record.upstreamEnd||f.c!==record.chrom||f.st!==record.strand||!/^[ACGTRYSWKMBDHVN]*$/.test(f.l+f.r))throw new Error('启动子外侧数据与所选转录本坐标不一致，未生成引物。');
  const sequence=f.l+record.upstream+f.r,start=f.l.length+record.upstream.length-target.length;
  return {...base,sequence,start,end:start+target.length,source:'同版本参考基因组的连续片段（已统一到基因转录方向）'};
 }
 return {...base,reason:'mRNA 模式仅设计所提供转录本的精确边界引物；CDS 两侧备用方案请切换为 CDS。'};
}

function countExact(sequence,query){let n=0,p=sequence.indexOf(query);while(p>=0){n++;p=sequence.indexOf(query,p+1);}return n;}
function annotatePair(pair,context,kind){
 const {start,end,sequence}=context;
 const f=pair.forward,r=pair.reverse;
 const product=sequence.slice(f.start,r.end);
 for(const [p,reverse] of [[f,false],[r,true]]){
  p.relativeStart=p.start-start+1;p.relativeEnd=p.end-start;
  // No coordinate zero: a base immediately before target position 1 is −1.
  if(p.relativeStart<=0)p.relativeStart--;if(p.relativeEnd<=0)p.relativeEnd--;
  p.exactSites=countExact(sequence,reverse?reverseComplement(p.sequence):p.sequence);
 }
 return {...pair,scheme:kind,productSequence:product,productSize:product.length,productWithTails:product.length+f.tail.length+r.tail.length,productWithTailsSequence:f.tail+product+reverseComplement(r.tail),leftExtension:start-f.start,rightExtension:r.end-end};
}
function flankCandidates(context,side,options,window){
 const {sequence,start,end}=context,min=options.minLen??18,max=options.maxLen??30,target=options.target??60;
 const tail=(side===0?options.forwardTail:options.reverseTail||'')?.replace(/\s/g,'').toUpperCase()||'';
 const lo=side===0?Math.max(0,start-window):end,hi=side===0?start:Math.min(sequence.length,end+window),out=[];
 for(let p=lo;p<=hi-min;p++)for(let len=min;len<=max&&p+len<=hi;len++){
  const seq=side===0?sequence.slice(p,p+len):reverseComplement(sequence.slice(p,p+len));
  if(!/^[ACGT]+$/.test(seq)||hasTandemRepeat(seq))continue;
  const tm=meltingTemp(seq,options);if(Math.abs(tm-target)>5)continue;
  const gc=100*(seq.match(/[GC]/g)||[]).length/len,homopolymer=Math.max(...seq.match(/(.)\1*/g).map(x=>x.length));
  if(gc<25||gc>75||homopolymer>5)continue;
  const fullSequence=tail+seq,self=complementarity(fullSequence,fullSequence);if(self.end>=5)continue;
  const hairpin=hairpinStem(fullSequence);if(hairpin>=6)continue;
  const distance=side===0?start-(p+len):p-end;
  out.push({sequence:seq,fullSequence,tail,length:len,tm,gc,homopolymer,self,hairpin,start:p,end:p+len,rank:Math.abs(tm-target)+Math.abs(len-22)*.15+self.end+distance*.02});
 }
 return out.sort((a,b)=>a.rank-b.rank||a.start-b.start||a.length-b.length);
}
export function designCloningSchemes(record,kind='cds',length=2000,options={},override=''){
 const context=cloningContext(record,kind,length,override),window=options.window??200;
 if(!Number.isInteger(window)||window<35||window>200)throw new Error('外侧搜索范围请填 35–200 bp（每侧）。');
 // This validates shared parameters even when no candidate can pass the filters.
 const fixed=designPrimers(context.target,{...options,avoidTandem:true});
 const primary={key:'exact',title:'首选方案 · 精确边界',pairs:fixed.pairs.map(p=>annotatePair({...p,forward:{...p.forward,start:context.start,end:context.start+p.forward.length},reverse:{...p.reverse,start:context.end-p.reverse.length,end:context.end}},context,'exact')),reason:fixed.reason};
 const backup={key:'outer',title:'备用方案 · 两侧外扩',pairs:[],reason:context.reason};
 const availableLeft=context.start,availableRight=context.sequence.length-context.end,min=options.minLen??18;
 if(!backup.reason&&(availableLeft<min||availableRight<min))backup.reason=`外侧序列不足：5′侧 ${availableLeft} bp，3′侧 ${availableRight} bp；每侧至少需 ${min} bp 才能把引物完全放在目标外。可在下方补充同一目标的两侧序列。`;
 if(!backup.reason){
  const forward=flankCandidates(context,0,options,window),reverse=flankCandidates(context,1,options,window),best=[];
  for(const f of forward){
   if(best.length===3&&f.rank+(reverse[0]?.rank??Infinity)>=best[2].score)break;
   for(const r of reverse){
    const lower=f.rank+r.rank;if(best.length===3&&lower>=best[2].score)break;
    const delta=Math.abs(f.tm-r.tm);if(delta>3)continue;
    if(best.length===3&&lower+delta>=best[2].score)continue;
    const cross=complementarity(f.fullSequence,r.fullSequence);if(cross.end>=5)continue;
    const score=lower+delta+cross.end;
    best.push({forward:{...f},reverse:{...r},tmDifference:delta,heterodimer:cross,score});best.sort((a,b)=>a.score-b.score);if(best.length>3)best.pop();
   }
  }
  backup.pairs=best.map(p=>annotatePair(p,context,'outer'));
  if(!best.length)backup.reason=`在目标两侧各最多 ${window} bp 内，未找到同时通过长度、Tm、GC 及结构初筛的外侧引物；未缩短完整目标。可调整参数或核对两侧序列。`;
 }
 return {schemes:[primary,backup],parameters:{...fixed.parameters,window},templateAmbiguous:fixed.templateAmbiguous,targetLength:context.target.length,kind,contextSource:context.source,availableLeft,availableRight,notes:kind==='promoter'?['启动子定义为注释 CDS 起点前的目标区，不等同于已验证的 TSS 启动子。备用产物会带入目标两侧的额外基因组序列，可能包含 CDS。']:kind==='cds'?['两套方案均保留完整 ATG 与末端终止密码子。备用产物额外含两侧 UTR，使用 cDNA 扩增。']:['mRNA 模式使用 cDNA 扩增。']};
}

export function cloningOrderTSV(result,record,firstOnly=true){
 const rows=[['名称','方案','序列_5to3','退火区','Tm_C','GC_percent','退火位置_目标相对','产物_bp','含附加序列_bp','参考版本','基因编号','转录本','验证状态']];
 for(const scheme of result.schemes)for(const [i,pair] of (firstOnly?scheme.pairs.slice(0,1):scheme.pairs).entries())for(const [d,p] of [['F',pair.forward],['R',pair.reverse]])rows.push([`${record.id}_${result.kind}_${scheme.key}_${i+1}_${d}`,scheme.title,p.fullSequence,p.sequence,p.tm.toFixed(2),p.gc.toFixed(1),`${p.relativeStart}..${p.relativeEnd}`,pair.productSize,pair.productWithTails,record.version,record.id,record.transcript,'计算初筛；未做全基因组特异性验证']);
 return '\uFEFF'+rows.map(r=>r.map(x=>String(x).replace(/[\t\r\n]/g,' ')).join('\t')).join('\n');
}
