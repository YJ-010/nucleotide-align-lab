import {reverseComplement,meltingTemp,complementarity,hairpinStem} from './clone-core.mjs?v=20261022-vector1';

const stopPattern=/T(?:AA|AG|GA)$/;
const IUPAC=/^[ACGTRYSWKMBDHVN]+$/;

export const VECTOR_TEMPLATES=Object.freeze([
 {id:'c-luc',name:'C-LUC（2b CDS）',short:'C-LUC',category:'荧光素酶融合',input:'cds',method:'同源重组',linearization:'KpnI / BamHI',sites:[['KpnI','GGTACC'],['BamHI','GGATCC']],forwardTail:'CGTCCCGGGGCGGTACC',reverseTail:'GTAGTCCATTTGTTGGATCCC',removeATG:false,removeStop:true,note:'保留起始密码子，去除终止密码子，用于 C 端 LUC 融合。'},
 {id:'n-luc',name:'N-LUC',short:'N-LUC',category:'荧光素酶融合',input:'cds',method:'同源重组',linearization:'KpnI / SalI',sites:[['KpnI','GGTACC'],['SalI','GTCGAC']],forwardTail:'GGGGGACGAGCTCGGTACC',reverseTail:'GTACGAGATCTGGTCGAC',removeATG:false,removeStop:true,note:'保留起始密码子，去除终止密码子；按模板连接 N-LUC。'},
 {id:'pri101-egfp',name:'PRI101-AN-EGFP（自定义）',short:'PRI101-EGFP',category:'绿色荧光融合',input:'cds',method:'同源重组',linearization:'NdeI / BamHI',sites:[['NdeI','CATATG'],['BamHI','GGATCC']],forwardTail:'CACTGTTGATACATATG',reverseTail:'CCCTTGCTCACGGATCC',removeATG:true,removeStop:true,note:'必须去除目标基因 ATG 和终止密码子；正向同源臂末端提供 ATG。'},
 {id:'1305',name:'1305（按模板）',short:'1305',category:'植物表达',input:'cds',method:'同源重组',linearization:'XbaI / BamHI',sites:[['XbaI','TCTAGA'],['BamHI','GGATCC']],forwardTail:'AGTCCGGAGCTAGCTCTAGA',reverseTail:'CCCTTGCTCACCATGGATCC',removeATG:false,removeStop:true,note:'按 Word 示例保留 ATG、去除终止密码子；该条由示例推定，订购前请核对实际 1305 载体图谱。',inferred:true},
 {id:'pgreen62',name:'pGreenII 62-SK',short:'pGreenII 62-SK',category:'植物表达',input:'cds',method:'同源重组',linearization:'XbaI / BamHI',sites:[['XbaI','TCTAGA'],['BamHI','GGATCC']],forwardTail:'GGTGGCGGCCGCTCTAGA',reverseTail:'TGCAGCCCGGGGGATCC',removeATG:false,removeStop:false,note:'按模板保留 ATG 和终止密码子。'},
 {id:'pgreen0800',name:'pGreen-0800（启动子）',short:'pGreen-0800',category:'启动子报告载体',input:'promoter',method:'同源重组',linearization:'XhoI / BamHI',sites:[['XhoI','CTCGAG'],['BamHI','GGATCC']],forwardTail:'CCGGGCCCCCCCTCGAG',reverseTail:'GCTCTAGAACTAGTGGATCC',removeATG:false,removeStop:false,note:'输入启动子 5′→3′ 序列，不做 ATG 或终止密码子处理。'},
 {id:'pgreen-gal4',name:'pGreenII 62-SK-Gal4',short:'62-SK-Gal4',category:'转录活性分析',input:'cds',method:'同源重组',linearization:'SmaI',sites:[['SmaI','CCCGGG']],forwardTail:'CTAGAACTAGTGGATCCCCCC',reverseTail:'ATATCGAATTCCTGCAGCCC',removeATG:false,removeStop:false,note:'按模板保留 ATG 和终止密码子；正向同源臂含额外 1 个 C 用于维持阅读框。'},
 {id:'pgreen-vp16',name:'pGreenII 62-SK-VP16',short:'62-SK-VP16',category:'转录抑制分析',input:'cds',method:'同源重组',linearization:'XbaI',sites:[['XbaI','TCTAGA']],forwardTail:'CCACCGCGGTGGCGGCCGCT',reverseTail:'GCCCGGGGGATCCACTAGTT',removeATG:false,removeStop:false,note:'按模板保留 ATG 和终止密码子；融合阅读框须结合实际载体图谱复核。'},
 {id:'pgadt7',name:'pGADT7',short:'pGADT7',category:'酵母双杂交',input:'cds',method:'同源重组',linearization:'EcoRI / BamHI',sites:[['EcoRI','GAATTC'],['BamHI','GGATCC']],forwardTail:'GCCATGGAGGCCAGTGAATTC',reverseTail:'GCAGCTCGAGCTCGATGGATCC',removeATG:false,removeStop:false,note:'保留终止密码子；表达框为 SV40 NLS–GAL4 AD–HA–目标蛋白，需核对阅读框。'}
]);

export function templateById(id){
 const value=VECTOR_TEMPLATES.find(x=>x.id===id);
 if(!value)throw new Error('请选择载体模板。');
 return value;
}

export function cleanDNA(raw,label='序列'){
 if(typeof raw!=='string'||!raw.trim())throw new Error('请粘贴'+label+'。');
 const lines=raw.replace(/\r/g,'').split('\n');
 const headers=lines.filter(line=>line.trim().startsWith('>'));
 if(headers.length>1)throw new Error(label+'输入框一次只能放一条序列。');
 const body=lines.filter(line=>!line.trim().startsWith('>')).join('');
 const sequence=body.toUpperCase().replace(/U/g,'T').replace(/[\s0-9]/g,'');
 if(!sequence)throw new Error('请粘贴'+label+'。');
 const invalid=[...new Set([...sequence].filter(x=>!'ACGTRYSWKMBDHVN'.includes(x)))];
 if(invalid.length)throw new Error(label+'含非法字符：'+invalid.join(' ')+'。请使用 DNA / IUPAC 碱基。');
 if(!IUPAC.test(sequence))throw new Error(label+'格式无法识别。');
 return sequence;
}

function cleanTail(value,label){
 const tail=String(value??'').toUpperCase().replace(/\s/g,'');
 if(!/^[ACGT]{8,60}$/.test(tail))throw new Error(label+'须为 8–60 个 A/C/G/T 碱基。');
 return tail;
}

export function prepareConstruct({templateId,cds='',promoter='',overrides={}}){
 const template=templateById(templateId);
 const source=cleanDNA(template.input==='promoter'?promoter:cds,template.input==='promoter'?'启动子序列':'CDS 序列');
 if(source.length<60||source.length>12000)throw new Error('目标序列须为 60–12,000 bp，当前为 '+source.length+' bp。');
 const removeATG=overrides.removeATG??template.removeATG;
 const removeStop=overrides.removeStop??template.removeStop;
 const forwardTail=cleanTail(overrides.forwardTail??template.forwardTail,'正向同源臂');
 const reverseTail=cleanTail(overrides.reverseTail??template.reverseTail,'反向同源臂');
 const warnings=[];
 let insert=source,removedATG=false,removedStop=false;
 if(template.input==='cds'){
  const hasATG=source.startsWith('ATG'),hasStop=stopPattern.test(source);
  if(source.length%3)warnings.push('原始 CDS 长度 '+source.length+' bp 不是 3 的倍数，请先核对开放阅读框。');
  if(removeATG){
   if(!hasATG)throw new Error('该载体规则要求去除起始密码子，但输入 CDS 不以 ATG 开头。请核对序列或在高级设置中更改。');
   insert=insert.slice(3);removedATG=true;
  }else if(!hasATG)warnings.push('输入 CDS 不以 ATG 开头；当前模板会原样保留 5′端。');
  if(removeStop){
   if(!hasStop)throw new Error('该载体规则要求去除终止密码子，但输入 CDS 不以 TAA/TAG/TGA 结尾。请核对序列或在高级设置中更改。');
   insert=insert.slice(0,-3);removedStop=true;
  }else if(!hasStop)warnings.push('输入 CDS 末端没有标准终止密码子；当前模板会原样保留 3′端。');
  if(insert.length%3)warnings.push('处理后插入片段为 '+insert.length+' bp，不是 3 的倍数，融合表达前必须复核阅读框。');
 }
 const ambiguous=(insert.match(/[^ACGT]/g)||[]).length;
 if(ambiguous)warnings.push('插入片段含 '+ambiguous+' 个不确定碱基；若位于两端会影响自动选引物。');
 const internalSites=[];
 for(const [name,site] of template.sites){
  let count=0,from=0,index;
  while((index=insert.indexOf(site,from))!==-1){count++;from=index+1;}
  if(count)internalSites.push({name,site,count});
 }
 if(internalSites.length)warnings.push('插入片段含线性化相关位点：'+internalSites.map(x=>x.name+' '+x.site+'×'+x.count).join('、')+'。同源重组通常不受影响，但若还要酶切 PCR 产物必须复核。');
 if(template.inferred)warnings.push('1305 的起止密码子规则是根据 Word 示例推定，必须与您手中的载体图谱再次核对。');
 return {template,source,insert,forwardTail,reverseTail,removeATG,removeStop,removedATG,removedStop,ambiguous,internalSites,warnings};
}

function longestRun(seq){return Math.max(...(seq.match(/(.)\1*/g)||['']).map(x=>x.length));}
function primerMetric(sequence,tail,options){
 const gc=100*[...sequence].filter(x=>'GC'.includes(x)).length/sequence.length;
 const fullSequence=tail+sequence;
 return {sequence,tail,fullSequence,length:sequence.length,totalLength:fullSequence.length,tm:meltingTemp(sequence,options),gc,homopolymer:longestRun(sequence),hairpin:hairpinStem(fullSequence),self:complementarity(fullSequence,fullSequence)};
}
function sideCandidates(insert,side,tail,options){
 const values=[];
 for(let length=options.minLen;length<=options.maxLen;length++){
  const sequence=side==='forward'?insert.slice(0,length):reverseComplement(insert.slice(-length));
  if(!/^[ACGT]+$/.test(sequence))continue;
  const metric=primerMetric(sequence,tail,options);
  metric.score=Math.abs(metric.tm-options.target)*3+Math.abs(metric.gc-50)/8+Math.abs(metric.length-22)*.15+Math.max(0,metric.homopolymer-4)*4+Math.max(0,metric.hairpin-4)*4+Math.max(0,metric.self.end-3)*6;
  values.push(metric);
 }
 return values.sort((a,b)=>a.score-b.score);
}
function pairWarnings(pair,target){
 const warnings=[];
 for(const [label,p] of [['正向',pair.forward],['反向',pair.reverse]]){
  if(Math.abs(p.tm-target)>5)warnings.push(label+'退火区 Tm 与目标相差 '+Math.abs(p.tm-target).toFixed(1)+' °C');
  if(p.gc<25||p.gc>75)warnings.push(label+'退火区 GC 为 '+p.gc.toFixed(1)+'%');
  if(p.homopolymer>5)warnings.push(label+'退火区含 '+p.homopolymer+' 个连续相同碱基');
  if(p.hairpin>=6)warnings.push(label+'完整引物简化发夹茎为 '+p.hairpin+' bp');
  if(p.self.end>=5)warnings.push(label+'完整引物 3′自互补为 '+p.self.end+' bp');
 }
 if(pair.tmDifference>3)warnings.push('两条引物 ΔTm 为 '+pair.tmDifference.toFixed(1)+' °C');
 if(pair.heterodimer.end>=5)warnings.push('两条完整引物 3′连续互补为 '+pair.heterodimer.end+' bp');
 return warnings;
}

export function designVectorPrimers({templateId,geneName='目标基因',cds='',promoter='',overrides={},options={}}){
 const prepared=prepareConstruct({templateId,cds,promoter,overrides});
 const settings={target:Number(options.target??60),minLen:Number(options.minLen??18),maxLen:Number(options.maxLen??30),salt:Number(options.salt??50),mg:Number(options.mg??1.5),dntp:Number(options.dntp??.2),concentration:Number(options.concentration??250)};
 if(!Number.isFinite(settings.target)||settings.target<50||settings.target>70)throw new Error('目标 Tm 请设为 50–70 °C。');
 if(!Number.isInteger(settings.minLen)||!Number.isInteger(settings.maxLen)||settings.minLen<18||settings.maxLen>35||settings.minLen>settings.maxLen)throw new Error('退火区长度范围须为 18–35 nt。');
 const forward=sideCandidates(prepared.insert,'forward',prepared.forwardTail,settings);
 const reverse=sideCandidates(prepared.insert,'reverse',prepared.reverseTail,settings);
 if(!forward.length||!reverse.length)throw new Error('目标片段的一端含不确定碱基，18–35 nt 范围内无法生成可靠退火区。请先核对序列两端。');
 const pairs=[];
 for(const f of forward.slice(0,7))for(const r of reverse.slice(0,7)){
  const tmDifference=Math.abs(f.tm-r.tm),heterodimer=complementarity(f.fullSequence,r.fullSequence);
  const pair={forward:f,reverse:r,tmDifference,heterodimer,productSize:prepared.insert.length,productWithTails:prepared.insert.length+prepared.forwardTail.length+prepared.reverseTail.length};
  pair.warnings=pairWarnings(pair,settings.target);
  pair.score=f.score+r.score+tmDifference*3+Math.max(0,heterodimer.end-3)*7+pair.warnings.length*2;
  pairs.push(pair);
 }
 pairs.sort((a,b)=>a.score-b.score||a.forward.length-b.forward.length||a.reverse.length-b.reverse.length);
 const safeName=(geneName.trim()||'目标基因').replace(/[\s>|/\\]+/g,'_').slice(0,80);
 const amplicon=prepared.forwardTail+prepared.insert+reverseComplement(prepared.reverseTail);
 return {...prepared,geneName:safeName,settings,pairs:pairs.slice(0,3),amplicon};
}

