import test from 'node:test';
import assert from 'node:assert/strict';
import {promoterSequence,cloningTemplate,compareDNA,meltingTemp,complementarity,hairpinStem,designPrimers,bindingCheck,reverseComplement,shardFor} from '../clone-core.mjs';
test('upstream orientation, plus/minus coordinates, chromosome-edge clipping and 2kb limit',()=>{
 const seq='A'.repeat(1000)+'C'.repeat(1000),plus={upstream:seq,strand:'+',upstreamStart:101,upstreamEnd:2100},minus={...plus,strand:'-'};
 assert.deepEqual(promoterSequence(plus,1000),{sequence:'C'.repeat(1000),start:1101,end:2100,requested:1000,truncated:false});
 assert.deepEqual(promoterSequence(minus,1000),{sequence:'C'.repeat(1000),start:101,end:1100,requested:1000,truncated:false});
 const short=promoterSequence({...minus,upstream:'ACGT'.repeat(30),upstreamEnd:220},2000);assert.equal(short.sequence.length,120);assert.equal(short.truncated,true);assert.equal(short.end,220);
 assert.throws(()=>promoterSequence(plus,2001),/60–2000/);assert.throws(()=>promoterSequence(plus,1.5),/60–2000/);
});
test('stop removal is explicit and preserves ATG',()=>{
 assert.equal(cloningTemplate({cds:'ATGCCCTAA'},'cds',2000,true).sequence,'ATGCCC');assert.equal(cloningTemplate({cds:'ATGCCCTAA'},'cds').sequence,'ATGCCCTAA');assert.throws(()=>cloningTemplate({cds:'ATGCCC'},'cds',2000,true),/终止/);
});
test('identity counts gaps and ambiguity; source sequences roundtrip through alignment',()=>{
 let r=compareDNA('ACGTN','ACGTN');assert.equal(r.identity,80);assert.equal(r.ambiguous,1);
 r=compareDNA('ACGT','ACGTT');assert.equal(r.identity,80);assert.equal(r.gaps,1);assert.equal(r.alignedA.replace(/-/g,''),'ACGT');assert.equal(r.alignedB.replace(/-/g,''),'ACGTT');
 r=compareDNA('GGGGACGTCCCC','ACGT');assert.ok(r.identity<40);assert.equal(r.coverageA,100/3);assert.equal(r.coverageB,100);assert.throws(()=>compareDNA('A'.repeat(7000),'A'.repeat(7000)),/上限/);
});
test('Tm matches independent Biopython DNA_NN3/saltcorr=5 fixtures, including symmetry',()=>{
 for(const [s,t] of [['CGTTCCAAAGATGTGGGCATGAGCTTAC',68.65341974660589],['ATGACCGTACAGGATCTGCA',61.50445853325266],['GCGCATATGCGC',50.91807925354436]])assert.ok(Math.abs(meltingTemp(s)-t)<1e-9);
 assert.throws(()=>meltingTemp('AANN'),/明确/);assert.throws(()=>meltingTemp('ACGT',{concentration:0}),/无效/);
});
test('antiparallel dimers anchor either 3-prime end; hairpins require an unpaired loop',()=>{
 assert.equal(complementarity('AAAAAA','TTTTTT').end,6);assert.equal(complementarity('AAAAAA','AAAAAA').any,0);assert.ok(complementarity('GGGGAAAAAA','TTTTTTCCCC').end>=6);assert.ok(hairpinStem('GCGCGCAAAAAAGCGCGC')>=6);assert.equal(hairpinStem('AAAAAAAAAAAA'),0);
});
test('cloning primer endpoints, reverse strand, tail lengths and refusal of low complexity',()=>{
 const seq='ATGACCGTACAGGATCTGCA'+'ACGT'.repeat(60)+'TGCATCAGGTACCGATGCTA',result=designPrimers(seq,{target:60});assert.ok(result.pairs.length>0);
 for(const p of result.pairs){assert.equal(p.forward.sequence,seq.slice(0,p.forward.length));assert.equal(p.reverse.sequence,reverseComplement(seq.slice(-p.reverse.length)));assert.equal(p.productSize,seq.length);assert.ok(p.tmDifference<=3);assert.equal(bindingCheck(p.forward.sequence,seq).mismatches,0);assert.equal(bindingCheck(p.reverse.sequence,seq,true).mismatches,0);}
 const tailed=designPrimers(seq,{forwardTail:'GAATTC',reverseTail:'TCTAGA'});for(const p of tailed.pairs){assert.ok(p.forward.fullSequence.startsWith('GAATTC'));assert.equal(p.productWithTails,seq.length+12);}
 assert.equal(designPrimers('N'.repeat(30)+'ACGT'.repeat(60)).pairs.length,0);assert.equal(designPrimers('A'.repeat(100)).pairs.length,0);assert.throws(()=>designPrimers(seq,{forwardTail:'G$TC'}),/附加/);
});
test('binding reports mismatch counts, 3-prime differences, coordinates and ties',()=>{
 const primer='AGTCCGTATG';let r=bindingCheck(primer,'AAA'+primer+'CCC');assert.equal(r.start,4);assert.equal(r.end,13);assert.equal(r.mismatches,0);r=bindingCheck(primer,'AAAAGTCCGTATAAAA');assert.equal(r.endMismatches,1);r=bindingCheck(primer,primer+'AAA'+primer);assert.equal(r.matches,2);assert.equal(bindingCheck(primer,'ACGT'),null);assert.equal(shardFor('Actinidia03760.t1'),shardFor('ACTINIDIA03760.T1'));
});
test('manual inputs retain unknown coordinates and missing sequence errors',()=>{
 const p=promoterSequence({upstream:'ACGT'.repeat(200),strand:'+',upstreamStart:null,upstreamEnd:null},500);assert.equal(p.start,null);assert.equal(p.end,null);assert.equal(p.sequence.length,500);
 assert.throws(()=>cloningTemplate({cds:''},'cds'),/未提供/);assert.throws(()=>cloningTemplate({upstream:''},'promoter'),/未提供/);
});
