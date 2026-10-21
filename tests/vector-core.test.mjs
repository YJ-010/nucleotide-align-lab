import test from 'node:test';
import assert from 'node:assert/strict';
import {VECTOR_TEMPLATES,cleanDNA,prepareConstruct,designVectorPrimers} from '../vector-core.mjs';
import {reverseComplement} from '../clone-core.mjs';

const cds='ATG'+Array.from({length:30},(_,i)=>['GCT','GAA','CAG','ATC','GAC','GGT','CTG','AAC'][i%8]).join('')+'TAA';
const promoter='ACGT'.repeat(40);

test('all Word templates have usable DNA homology arms',()=>{
 assert.equal(VECTOR_TEMPLATES.length,9);
 for(const item of VECTOR_TEMPLATES){
  assert.match(item.forwardTail,/^[ACGT]{8,60}$/);
  assert.match(item.reverseTail,/^[ACGT]{8,60}$/);
  assert.ok(['cds','promoter'].includes(item.input));
 }
});

test('PRI101 removes ATG and terminal stop exactly once',()=>{
 const result=prepareConstruct({templateId:'pri101-egfp',cds,promoter});
 assert.equal(result.insert,cds.slice(3,-3));
 assert.equal(result.removedATG,true);
 assert.equal(result.removedStop,true);
});

test('pGreenII 62-SK retains ATG and stop',()=>{
 const result=prepareConstruct({templateId:'pgreen62',cds});
 assert.equal(result.insert,cds);
 assert.equal(result.removedATG,false);
 assert.equal(result.removedStop,false);
});

test('pGreen-0800 uses promoter rather than CDS',()=>{
 const result=prepareConstruct({templateId:'pgreen0800',cds:'ATG'+'GCT'.repeat(30)+'TAA',promoter});
 assert.equal(result.source,promoter);
 assert.equal(result.insert,promoter);
});

test('primer designer joins template tails in 5-prime direction',()=>{
 const result=designVectorPrimers({templateId:'c-luc',geneName:'AcERF1',cds});
 const pair=result.pairs[0];
 assert.ok(pair.forward.fullSequence.startsWith(result.forwardTail));
 assert.ok(pair.reverse.fullSequence.startsWith(result.reverseTail));
 assert.ok(result.amplicon.startsWith(result.forwardTail));
 assert.ok(result.amplicon.endsWith(reverseComplement(result.reverseTail)));
 assert.equal(result.insert,cds.slice(0,-3));
 assert.equal(result.pairs.length,3);
});

test('related internal restriction site is reported',()=>{
 const withSite='ATG'+'GCT'.repeat(10)+'GGATCC'+'GAA'.repeat(18)+'TAA';
 const result=prepareConstruct({templateId:'c-luc',cds:withSite});
 assert.ok(result.internalSites.some(x=>x.name==='BamHI'));
});

test('multiple FASTA records are rejected instead of concatenated',()=>{
 assert.throws(()=>cleanDNA('>a\nACGT\n>b\nACGT','CDS'),/一次只能放一条/);
});

