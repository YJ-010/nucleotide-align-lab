import test from 'node:test';
import assert from 'node:assert/strict';
import {cloningContext,designCloningSchemes,cloningOrderTSV} from '../cloning-design.mjs';
import {reverseComplement} from '../clone-core.mjs';

const cds='ATGACCGTACAGGATCTGCA'+'ACGT'.repeat(60)+'CC'+'TGCATCAGGTACCGATGCTATAA';
let seed=47;function dna(n){let s='';for(let i=0;i<n;i++){seed=(Math.imul(seed,1664525)+1013904223)>>>0;s+='ACGT'[seed>>>30];}return s;}
const left=dna(200),right=dna(200);
const record={id:'PLANT',transcript:'PLANT.t1',version:'hy3',cds,mrna:left+cds+right,mrnaInfo:{cdsStart:201,cdsEnd:200+cds.length}};
test('both schemes retain complete CDS and use correct reverse primer orientation',()=>{
 const r=designCloningSchemes(record);assert.equal(cds.length%3,0);
 assert.ok(r.schemes.every(s=>s.pairs.length>0));
 for(const scheme of r.schemes)for(const p of scheme.pairs){
  assert.equal(p.productSequence,record.mrna.slice(p.forward.start,p.reverse.end));
  assert.ok(p.productSequence.includes(cds));
  assert.equal(p.reverse.sequence,reverseComplement(record.mrna.slice(p.reverse.start,p.reverse.end)));
  assert.equal(p.productSize,cds.length+p.leftExtension+p.rightExtension);
  if(scheme.key==='exact'){
   assert.ok(p.forward.sequence.startsWith('ATG'));assert.ok(p.reverse.sequence.startsWith('TTA'));
   assert.equal(p.productSequence,cds);assert.equal(p.forward.relativeStart,1);assert.equal(p.reverse.relativeEnd,cds.length);
  }else{
   assert.ok(p.forward.end<=200);assert.ok(p.reverse.start>=200+cds.length);
   assert.ok(p.forward.relativeEnd<0);assert.ok(p.reverse.relativeStart>cds.length);
   assert.ok(p.leftExtension<=200&&p.rightExtension<=200);
  }
 }
});
test('UTR absence never substitutes genomic sequence; one missing flank blocks only backup',()=>{
 const bare={...record,mrna:'',mrnaInfo:null,upstream:left};const r=designCloningSchemes(bare);
 assert.ok(r.schemes[0].pairs.length);assert.equal(r.schemes[1].pairs.length,0);assert.match(r.schemes[1].reason,/UTR/);
 const one=designCloningSchemes({...record,mrna:left+cds});assert.equal(one.schemes[1].pairs.length,0);assert.match(one.schemes[1].reason,/3′侧 0 bp/);
});
test('precise boundary failure still allows valid outer primers without target trimming',()=>{
 const difficult='ATG'+'A'.repeat(300)+'TAA';const r=designCloningSchemes({...record,cds:difficult,mrna:left+difficult+right,mrnaInfo:{cdsStart:201,cdsEnd:200+difficult.length}});
 assert.equal(r.schemes[0].pairs.length,0);assert.ok(r.schemes[1].pairs.length);
 assert.ok(r.schemes[1].pairs.every(p=>p.productSequence.includes(difficult)));
});
test('promoter primary is 2000 bp and backup uses contiguous genomic flanks on either strand',()=>{
 const target=cds+dna(2000-cds.length);
 for(const strand of ['+','-']){
  const r={...record,strand,chrom:'chr1',upstream:target,upstreamStart:501,upstreamEnd:2500,promoterFlanks:{l:left,r:right,s:501,e:2500,n:2000,c:'chr1',st:strand}};
  const c=cloningContext(r,'promoter');assert.equal(c.sequence,left+target+right);
  const result=designCloningSchemes(r,'promoter');assert.equal(result.targetLength,2000);assert.ok(result.schemes[1].pairs.length);
  for(const p of result.schemes[0].pairs)assert.equal(p.productSequence,target);
  for(const p of result.schemes[1].pairs){assert.ok(p.productSequence.includes(target));assert.equal(p.productSequence,(left+target+right).slice(p.forward.start,p.reverse.end));}
  const short=cloningContext(r,'promoter',1000);assert.equal(short.target,target.slice(-1000));assert.equal(short.start,1200);
  assert.throws(()=>cloningContext({...r,promoterFlanks:{...r.promoterFlanks,st:strand==='+'?'-':'+'}},'promoter'),/不一致/);
 }
 assert.throws(()=>cloningContext({upstream:left,strand:'+'},'promoter'),/完整 2000/);
});
test('manual context must uniquely and exactly contain full target; incomplete CDS rejected',()=>{
 assert.equal(cloningContext({...record,mrna:''},'cds',2000,left+cds+right).start,200);
 assert.throws(()=>cloningContext(record,'cds',2000,left+cds.slice(3)+right),/完整包含/);
 assert.throws(()=>cloningContext(record,'cds',2000,cds+cds),/仅一次/);
 assert.throws(()=>cloningContext({...record,cds:cds.slice(3)}),/ATG/);
 assert.throws(()=>cloningContext({...record,cds:cds.slice(0,-3)}),/终止/);
 assert.throws(()=>designCloningSchemes(record,'cds',2000,{window:201}),/35–200/);
});
test('tails are included with correct orientation in products; TSV exports one pair per scheme',()=>{
 const r=designCloningSchemes(record,'cds',2000,{forwardTail:'GAATTC',reverseTail:'AGATCT'});
 assert.ok(r.schemes.every(s=>s.pairs.length));
 for(const s of r.schemes)for(const p of s.pairs){assert.equal(p.productWithTailsSequence,'GAATTC'+p.productSequence+reverseComplement('AGATCT'));assert.equal(p.productWithTails,p.productSize+12);}
 const rows=cloningOrderTSV(r,record).trim().split('\n');assert.equal(rows.length,5);assert.match(rows[1],/PLANT_cds_exact_1_F/);assert.match(rows[3],/PLANT_cds_outer_1_F/);
});
