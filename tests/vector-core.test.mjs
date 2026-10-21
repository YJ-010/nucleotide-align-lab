import test from 'node:test';
import assert from 'node:assert/strict';
import {VECTOR_TEMPLATES,cleanDNA,prepareConstruct,designVectorPrimers} from '../vector-core.mjs';
import {reverseComplement,vectorInputFromRecord,cloningTemplate} from '../clone-core.mjs';
import {vectorDemo} from '../vector-demo.mjs';

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

test('clone-to-vector transfer keeps original CDS even after clone stop removal',()=>{
 const record={version:'hy2',id:'test',cds,upstream:promoter};
 const processed=cloningTemplate(record,'cds',2000,true);
 assert.equal(processed.sequence,cds.slice(0,-3));
 const input=vectorInputFromRecord(record,'cds',80);
 assert.equal(input.cds,cds);
 assert.equal(input.promoter,promoter.slice(-80));
 assert.equal(prepareConstruct({templateId:'pri101-egfp',...input}).insert,cds.slice(3,-3));
});

test('import does not substitute mRNA or reuse a missing promoter',()=>{
 const record={version:'white',id:'test',cds};
 assert.equal(vectorInputFromRecord(record,'cds').promoter,'');
 assert.throws(()=>vectorInputFromRecord(record,'mrna'),/mRNA/);
});

test('in-frame internal stop codons are warned without silently changing input',()=>{
 const internal='ATG'+'GCT'.repeat(10)+'TGA'+'GAA'.repeat(18)+'TAA';
 const result=prepareConstruct({templateId:'pgreen62',cds:internal});
 assert.equal(result.insert,internal);
 assert.ok(result.warnings.some(x=>x.includes('内部发现 1 个终止密码子')));
});

test('short inserts after codon removal are rejected',()=>{
 assert.throws(()=>prepareConstruct({templateId:'pri101-egfp',cds:'ATG'+'GCT'.repeat(18)+'TAA'}),/至少 60/);
});

test('all length combinations are evaluated and passing pairs are ranked first',()=>{
 const demo=vectorDemo();
 for(const templateId of ['pri101-egfp','c-luc','pgreen0800']){
  const result=designVectorPrimers({templateId,...demo});
  assert.equal(result.evaluatedPairs,13*13);
  assert.ok(result.passingPairs>0);
  assert.equal(result.pairs[0].warnings.length,0);
 }
});

test('fixed primer length returns one candidate rather than a promised three',()=>{
 const result=designVectorPrimers({templateId:'c-luc',...vectorDemo(),options:{minLen:22,maxLen:22}});
 assert.equal(result.evaluatedPairs,1);
 assert.equal(result.pairs.length,1);
});

test('repeated annealing sequences are flagged within the submitted template',()=>{
 const result=designVectorPrimers({templateId:'pgreen0800',promoter:'ACGT'.repeat(80)});
 assert.equal(result.passingPairs,0);
 assert.ok(result.pairs[0].warnings.some(x=>x.includes('完全匹配')));
});

test('all nine templates generate candidates using synthetic input',()=>{
 for(const template of VECTOR_TEMPLATES){
  const result=designVectorPrimers({templateId:template.id,...vectorDemo()});
  assert.ok(result.pairs.length>0);
  assert.equal(result.amplicon,result.forwardTail+result.insert+reverseComplement(result.reverseTail));
 }
});

