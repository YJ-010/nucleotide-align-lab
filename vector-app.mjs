import {VECTOR_TEMPLATES,templateById,designVectorPrimers} from './vector-core.mjs?v=20261021-kiwi2';
import {vectorDemo} from './vector-demo.mjs?v=20261021-kiwi2';

const $=id=>document.getElementById(id);
const element=(tag,text,className)=>{const node=document.createElement(tag);if(text!==undefined)node.textContent=text;if(className)node.className=className;return node;};
let currentResult=null,reportText='',orderText='';

function setStatus(text,error=false){
 const node=$('vector-status');
 node.textContent=text;
 node.classList.toggle('error',error);
 node.classList.toggle('clone-status-error',error);
}
function save(name,text,type='text/plain;charset=utf-8'){
 const link=document.createElement('a');
 link.href=URL.createObjectURL(new Blob([text],{type}));
 link.download=name;link.click();
 setTimeout(()=>URL.revokeObjectURL(link.href),1000);
}
function lines(sequence){return sequence.match(/.{1,80}/g)?.join('\n')||'';}
function inputLength(id){
 const value=$(id).value.replace(/^>.*$/gm,'').toUpperCase().replace(/U/g,'T').replace(/[\s0-9]/g,'');
 $(id+'-count').textContent=value.length?value.length+' bp':'尚未输入';
}

function buildTemplateSelect(){
 const groups=new Map();
 for(const item of VECTOR_TEMPLATES){
  if(!groups.has(item.category))groups.set(item.category,[]);
  groups.get(item.category).push(item);
 }
 for(const [category,items] of groups){
  const group=document.createElement('optgroup');group.label=category;
  for(const item of items){const option=element('option',item.name);option.value=item.id;group.append(option);}
  $('vector-template').append(group);
 }
}

function renderRule(reset=true){
 const item=templateById($('vector-template').value);
 if(reset){
  $('vector-forward-tail').value=item.forwardTail;
  $('vector-reverse-tail').value=item.reverseTail;
  $('vector-remove-atg').checked=item.removeATG;
  $('vector-remove-stop').checked=item.removeStop;
 }
 const source=item.input==='promoter'?'启动子序列':'CDS 序列';
 const removeATG=$('vector-remove-atg').checked,removeStop=$('vector-remove-stop').checked;
 const forwardTail=$('vector-forward-tail').value,reverseTail=$('vector-reverse-tail').value;
 const customized=removeATG!==item.removeATG||removeStop!==item.removeStop||forwardTail!==item.forwardTail||reverseTail!==item.reverseTail;
 $('vector-required-sequence').textContent='本载体读取：'+source;
 $('vector-cds-card').classList.toggle('required-sequence',item.input==='cds');
 $('vector-promoter-card').classList.toggle('required-sequence',item.input==='promoter');
 $('vector-remove-atg').disabled=item.input==='promoter';
 $('vector-remove-stop').disabled=item.input==='promoter';
 const box=$('vector-template-summary');box.replaceChildren();
 const title=element('div',undefined,'vector-rule-head');
 title.append(element('strong',item.name),element('span',item.method+' · '+item.linearization));
 box.append(title);
 const grid=element('dl',undefined,'vector-rule-grid');
 const rows=[
  ['输入对象',source],
  ['载体线性化',item.linearization],
  ['起始密码子',item.input==='promoter'?'不适用':removeATG?'去除 CDS 的 ATG':'保留'],
  ['终止密码子',item.input==='promoter'?'不适用':removeStop?'去除':'保留'],
  ['正向同源臂',forwardTail||'尚未填写'],
  ['反向同源臂',reverseTail||'尚未填写']
 ];
 for(const [term,value] of rows){grid.append(element('dt',term),element('dd',value));}
 box.append(grid,element('p',customized?'当前为自定义规则；上方展示的同源臂与密码子设置将用于计算。原模板说明：'+item.note:item.note,(item.inferred||customized)?'vector-rule-note warning':'vector-rule-note'));
 $('vector-design').textContent='按 '+item.short+' 模板生成引物 →';
 clearResults();
}

function clearResults(){
 currentResult=null;reportText='';orderText='';
 $('vector-results').hidden=true;
 $('vector-output').replaceChildren();
 $('vector-actions').hidden=true;
 setStatus('');
}

function primerSequence(primer){
 const code=element('code',undefined,'segmented-primer');
 code.append(element('span',primer.tail,'primer-tail'),element('span',primer.sequence,'primer-binding'));
 return code;
}
function statText(primer){
 return '总长 '+primer.totalLength+' nt · 退火区 '+primer.length+' nt · Tm '+primer.tm.toFixed(2)+' °C · GC '+primer.gc.toFixed(1)+'%';
}
function copyButton(sequence,label){
 const button=element('button','复制','copy-primer small-action');button.type='button';
 button.addEventListener('click',async()=>{
  try{await navigator.clipboard.writeText(sequence);button.textContent='已复制';setTimeout(()=>button.textContent='复制',1400);}
  catch{setStatus('浏览器未允许复制，请使用“下载订购清单”。',true);}
 });
 button.setAttribute('aria-label','复制'+label+'引物');
 return button;
}
function primerRow(label,primer){
 const row=element('div',undefined,'vector-primer-row');
 const meta=element('div',undefined,'vector-primer-label');meta.append(element('strong',label),element('small','5′→3′'));
 const value=element('div',undefined,'vector-primer-value');
 const top=element('div',undefined,'vector-primer-code');top.append(primerSequence(primer),copyButton(primer.fullSequence,label));
 value.append(top,element('p',statText(primer),'field-hint'));
 row.append(meta,value);return row;
}

function buildReport(result){
 const rows=[
  '猕猴桃基因工坊 · 载体构建与自动引物设计报告',
  '基因/片段：'+result.geneName,
  '载体：'+result.template.name,
  '构建方式：'+result.template.method,
  '线性化：'+result.template.linearization,
  '输入对象：'+(result.template.input==='promoter'?'启动子':'CDS'),
  '原始序列：'+result.source.length+' bp',
  '处理后插入片段：'+result.insert.length+' bp',
  'ATG：'+(result.template.input==='promoter'?'不适用':result.removedATG?'已去除':'保留'),
  '终止密码子：'+(result.template.input==='promoter'?'不适用':result.removedStop?'已去除':'保留'),
  '正向同源臂：'+result.forwardTail,
  '反向同源臂：'+result.reverseTail,
  '参数：目标 Tm '+result.settings.target+' °C，退火区 '+result.settings.minLen+'–'+result.settings.maxLen+' nt',
  'Tm 默认反应条件：单价盐 '+result.settings.salt+' mM，Mg²⁺ '+result.settings.mg+' mM，总 dNTP '+result.settings.dntp+' mM，寡核苷酸 '+result.settings.concentration+' nM',
  '筛查：比较 '+result.evaluatedPairs+' 组长度组合，'+result.passingPairs+' 组通过本工具引物初筛；非全基因组特异性验证。',
  '',
  '候选 1（5′→3′；'+(result.pairs[0].warnings.length?'未通过全部初筛，请复核后再订购':'通过本工具初筛，仍须核对实际载体')+'）',
  result.geneName+'_'+result.template.short+'_F：'+result.pairs[0].forward.fullSequence,
  result.geneName+'_'+result.template.short+'_R：'+result.pairs[0].reverse.fullSequence,
  'F：'+statText(result.pairs[0].forward),
  'R：'+statText(result.pairs[0].reverse),
  'ΔTm：'+result.pairs[0].tmDifference.toFixed(2)+' °C',
  ''
 ];
 if(result.warnings.length){rows.push('模板/序列提醒：',...result.warnings.map(x=>'· '+x),'');}
 if(result.pairs[0].warnings.length){rows.push('引物初筛提醒：',...result.pairs[0].warnings.map(x=>'· '+x),'');}
 rows.push('处理后插入片段 FASTA：','>'+result.geneName+'|'+result.template.id+'|processed_insert',lines(result.insert),'','含两端同源序列的 PCR 产物 FASTA：','>'+result.geneName+'|'+result.template.id+'|homology_extended_amplicon',lines(result.amplicon),'','说明：该结果按用户提供的 Word 模板自动拼接；订购前仍须结合实际载体图谱、阅读框和测序方向复核。');
 return rows.join('\n');
}
function buildOrder(result){
 const pair=result.pairs[0],base=result.geneName+'_'+result.template.short.replace(/\s+/g,'_');
 const rows=[['名称','序列_5to3','总长度_nt','退火区','退火区长度_nt','Tm_C','GC_percent','载体','线性化','处理规则','状态']];
 for(const [suffix,primer] of [['F',pair.forward],['R',pair.reverse]]){
  rows.push([base+'_'+suffix,primer.fullSequence,primer.totalLength,primer.sequence,primer.length,primer.tm.toFixed(2),primer.gc.toFixed(1),result.template.name,result.template.linearization,result.template.input==='promoter'?'起止密码子处理不适用':(result.removedATG?'去ATG；':'保留ATG；')+(result.removedStop?'去终止':'保留终止'),pair.warnings.length||result.warnings.length?'需复核；详见报告':'引物初筛通过；载体仍需复核']);
 }
 return '\uFEFF'+rows.map(row=>row.map(value=>{const s=String(value).replace(/[\t\r\n]/g,' ');return /^[=+@-]/.test(s)?"'"+s:s;}).join('\t')).join('\n');
}

function renderResults(result){
 const output=$('vector-output');output.replaceChildren();
 const summary=element('div',undefined,'vector-summary-grid');
 const facts=[
  ['模板',result.template.short],
  ['插入片段',result.insert.length+' bp'],
  ['候选 1 长度',result.pairs[0].forward.totalLength+' / '+result.pairs[0].reverse.totalLength+' nt'],
  ['初筛',result.pairs[0].warnings.length?'需复核 '+result.pairs[0].warnings.length+' 项':'通过']
 ];
 for(const [label,value] of facts){const card=element('div',undefined,'vector-summary-card');card.append(element('span',label),element('strong',value));summary.append(card);}
 output.append(summary);
 output.append(element('p','已比较 '+result.evaluatedPairs+' 组引物长度组合，'+result.passingPairs+' 组通过本工具初筛。展示排名靠前的 '+result.pairs.length+' 组；通过初筛不等于实验验证通过。','field-hint'));
 const decision=element('section',undefined,'vector-decision');
 decision.append(element('h3','自动处理结果'));
 const pills=element('div',undefined,'vector-decision-pills');
 const rules=[
  result.template.input==='promoter'?'使用启动子':'使用 CDS',
  result.template.input==='promoter'?'不处理 ATG':result.removedATG?'去除 ATG':'保留 ATG',
  result.template.input==='promoter'?'不处理终止密码子':result.removedStop?'去除终止密码子':'保留终止密码子',
  '线性化 '+result.template.linearization
 ];
 for(const value of rules)pills.append(element('span',value));
 decision.append(pills,element('p',result.template.note,'field-hint'));
 output.append(decision);
 if(result.warnings.length){
  const warnings=element('section',undefined,'vector-warnings');warnings.append(element('h3','构建前提醒'));
  const list=element('ul');for(const warning of result.warnings)list.append(element('li',warning));warnings.append(list);output.append(warnings);
 }
 const alternatives=element('details',undefined,'vector-alternatives');
 alternatives.append(element('summary','展开其他 '+Math.max(0,result.pairs.length-1)+' 组备选引物'));
 result.pairs.forEach((pair,index)=>{
  const card=element('article',undefined,'vector-pair'+(index===0?' preferred':''));
  const heading=element('div',undefined,'vector-pair-heading');
  heading.append(element('div',undefined,'vector-pair-number'),element('h3',index===0?'候选 1 · '+(pair.warnings.length?'复核后再订购':'通过引物初筛'):'备选 '+(index+1)));
  heading.firstChild.textContent='0'+(index+1);
  heading.append(element('span',pair.warnings.length?'需复核':'初筛通过',pair.warnings.length?'vector-badge caution':'vector-badge'));
  card.append(heading,primerRow('F 正向',pair.forward),primerRow('R 反向',pair.reverse));
  const legend=element('p',undefined,'primer-color-legend');
  legend.append(element('span','同源臂','legend-tail'),document.createTextNode(' + '),element('span','基因退火区','legend-binding'),document.createTextNode('；Tm 只按退火区计算。'));
  card.append(legend,element('p','ΔTm '+pair.tmDifference.toFixed(2)+' °C · 两引物 3′连续互补 '+pair.heterodimer.end+' bp · PCR 插入片段 '+pair.productSize+' bp','field-hint'));
  if(pair.warnings.length){const list=element('ul',undefined,'primer-warning-list');for(const warning of pair.warnings)list.append(element('li',warning));card.append(list);}
  (index===0?output:alternatives).append(card);
 });
 if(result.pairs.length>1)output.append(alternatives);
 const sequences=element('details',undefined,'vector-sequence-preview');
 sequences.append(element('summary','查看处理后的插入片段与含同源臂 PCR 产物'));
 sequences.append(element('h3','插入片段 · '+result.insert.length+' bp'),element('pre',lines(result.insert)));
 sequences.append(element('h3','含同源臂 PCR 产物 · '+result.amplicon.length+' bp'),element('pre',lines(result.amplicon)));
 output.append(sequences);
 const next=element('section',undefined,'vector-next-step');
 next.append(element('h3','下一步怎么做'));
 const list=element('ol');
 for(const text of [result.pairs[0].warnings.length?'当前候选未通过全部初筛，请先处理风险提示、调整参数或用专业软件复核，暂勿直接订购。':'核对实际载体图谱、阅读框与引物特异性，通过复核后再订购 F / R 引物。','按实验室 SOP 扩增目标片段、线性化载体并完成同源重组。','挑克隆测序；下方可把预期插入片段送入测序比对。该参照不含完整载体，不能单独验证两端连接位点。'])list.append(element('li',text));
 next.append(list);output.append(next);
 $('vector-results').hidden=false;$('vector-actions').hidden=false;
 reportText=buildReport(result);orderText=buildOrder(result);
 $('vector-results').scrollIntoView({behavior:'smooth',block:'start'});
}

function design(){
 clearResults();
 try{
  const result=designVectorPrimers({
   templateId:$('vector-template').value,
   geneName:$('vector-gene-name').value,
   cds:$('vector-cds').value,
   promoter:$('vector-promoter').value,
   overrides:{removeATG:$('vector-remove-atg').checked,removeStop:$('vector-remove-stop').checked,forwardTail:$('vector-forward-tail').value,reverseTail:$('vector-reverse-tail').value},
   options:{target:Number($('vector-tm').value),minLen:Number($('vector-min').value),maxLen:Number($('vector-max').value)}
  });
  currentResult=result;renderResults(result);
  setStatus('已按 '+result.template.name+' 生成 '+result.pairs.length+' 组候选。'+(result.passingPairs?'通过引物初筛的组合已优先排列。':'没有组合通过全部引物初筛，请先复核，勿直接订购。')+(result.warnings.length?' 请先查看构建前提醒。':''));
 }catch(error){setStatus(error.message,true);}
}

function loadExample(){
 const demo=vectorDemo();
 $('vector-gene-name').value=demo.geneName;
 $('vector-cds').value='>'+demo.geneName+'_CDS\n'+demo.cds;
 $('vector-promoter').value='>'+demo.geneName+'_promoter\n'+demo.promoter;
 $('vector-template').value='pri101-egfp';
 renderRule(true);inputLength('vector-cds');inputLength('vector-promoter');
 setStatus('已载入模拟示例，可直接点击生成。示例不代表真实猕猴桃基因。');
}
function clearAll(){
 for(const id of ['vector-gene-name','vector-cds','vector-promoter'])$(id).value='';
 inputLength('vector-cds');inputLength('vector-promoter');clearResults();
}

buildTemplateSelect();
$('vector-template').value='c-luc';
renderRule(true);
$('vector-template').addEventListener('change',()=>renderRule(true));
for(const id of ['vector-cds','vector-promoter'])$(id).addEventListener('input',()=>{inputLength(id);clearResults();});
for(const id of ['vector-gene-name','vector-tm','vector-min','vector-max'])$(id).addEventListener('input',clearResults);
for(const id of ['vector-forward-tail','vector-reverse-tail','vector-remove-atg','vector-remove-stop'])$(id).addEventListener('input',()=>renderRule(false));
$('vector-reset-rules').addEventListener('click',()=>{renderRule(true);setStatus('已恢复此载体的模板默认规则，输入序列未改变。');});
$('vector-design').addEventListener('click',design);
$('vector-example').addEventListener('click',loadExample);
$('vector-clear').addEventListener('click',clearAll);
$('vector-copy-all').addEventListener('click',async()=>{try{await navigator.clipboard.writeText(reportText);setStatus('引物与构建报告已复制。');}catch{setStatus('浏览器未允许复制，请点击下载报告。',true);}});
$('vector-download-report').addEventListener('click',()=>{if(reportText)save('KiwiGene_'+currentResult.geneName+'_'+currentResult.template.id+'_引物报告.txt',reportText);});
$('vector-download-tsv').addEventListener('click',()=>{if(orderText)save('KiwiGene_'+currentResult.geneName+'_'+currentResult.template.id+'_引物订购单.tsv',orderText,'text/tab-separated-values;charset=utf-8');});
$('vector-download-fasta').addEventListener('click',()=>{if(currentResult)save('KiwiGene_'+currentResult.geneName+'_'+currentResult.template.id+'_构建序列.fasta','>'+currentResult.geneName+'|processed_insert\n'+lines(currentResult.insert)+'\n>'+currentResult.geneName+'|homology_extended_amplicon\n'+lines(currentResult.amplicon)+'\n','text/plain;charset=utf-8');});
$('vector-send-to-alignment').addEventListener('click',()=>{
 if(!currentResult)return;
 window.dispatchEvent(new CustomEvent('seqlab-reference',{detail:{name:currentResult.geneName+'|'+currentResult.template.id+'|expected_insert',sequence:currentResult.insert}}));
 window.dispatchEvent(new CustomEvent('seqlab-show-workspace',{detail:{workspace:'sequence'}}));
 $('sequence-input').scrollIntoView({behavior:'smooth',block:'center'});
});
window.addEventListener('seqlab-vector-input',event=>{
 const detail=event.detail||{};
 $('vector-gene-name').value=detail.name||'目标基因';
 $('vector-cds').value=detail.cds??(detail.kind==='cds'?detail.sequence||'':'');
 $('vector-promoter').value=detail.promoter??(detail.kind==='promoter'?detail.sequence||'':'');
 if(detail.kind==='promoter')$('vector-template').value='pgreen0800';
 else if($('vector-template').value==='pgreen0800')$('vector-template').value='c-luc';
 renderRule(true);inputLength('vector-cds');inputLength('vector-promoter');
 window.dispatchEvent(new CustomEvent('seqlab-show-workspace',{detail:{workspace:'vector'}}));
 setStatus('已带入同一基因的'+($('vector-cds').value&&$('vector-promoter').value?'原始 CDS 和上游区':$('vector-cds').value?'原始 CDS':'上游区')+'。起止密码子统一由载体规则处理；缺少的序列不会沿用上一个基因。');
});

