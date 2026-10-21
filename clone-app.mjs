import {parseSequences} from './align-core.mjs?v=20261021-kiwi2';
import {VERSIONS,shardFor,cloningTemplate,compareDNA,designPrimers,bindingCheck,fastaText,withTranscript,sequenceDisplay,vectorInputFromRecord} from './clone-core.mjs?v=20261021-kiwi2';
const $=id=>document.getElementById(id), labels=Object.fromEntries(VERSIONS), DATA='./data/kiwi-20261019/';
let revision=0,manifest=null,source=null,selected=new Map(),candidateLists=new Map(),comparisons=[],activePair=null,primerResult=null,primerReport='',worker=null,jobId=0;
const jobs=new Map(),cache=new Map();
function el(tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n;}
function say(id,text,error=false){const n=$(id);n.textContent=text;n.classList.toggle('clone-status-error',error);}
function cancelJobs(){if(worker){worker.terminate();worker=null;}for(const job of jobs.values()){clearTimeout(job.timer);job.reject(new Error('参数已更改，请重新计算。'));}jobs.clear();}
function clearPrimer(){primerResult=null;primerReport='';$('primer-results').replaceChildren();$('primer-downloads').hidden=true;say('primer-status','');}
function invalidateComparison(){revision++;cancelJobs();comparisons=[];activePair=null;$('comparison-matrix').replaceChildren();$('pair-detail').hidden=true;say('compare-status','参数已更新，点击比较当前候选。');clearPrimer();}
function clearGenes(){invalidateComparison();source=null;selected=new Map();candidateLists=new Map();$('gene-results').hidden=true;$('gene-cards').replaceChildren();}
async function jsonFetch(path,gzip=false){
 const controller=new AbortController(),timeout=setTimeout(()=>controller.abort(),60000);
 try{const r=await fetch(DATA+path+(gzip?'':'?v=20261019-clone4'),{signal:controller.signal});if(!r.ok)throw new Error(`参考数据读取失败（HTTP ${r.status}），请稍后重试。`);if(!gzip)return await r.json();
 const b=new Uint8Array(await r.arrayBuffer());let decoded;
 if(b[0]===31&&b[1]===139){if(typeof DecompressionStream==='undefined')throw new Error('此浏览器不支持参考库解压，请使用近期版本的 Chrome、Edge、Firefox 或 Safari。');decoded=await new Response(new Blob([b]).stream().pipeThrough(new DecompressionStream('gzip'))).text();}else decoded=new TextDecoder().decode(b);
 return JSON.parse(decoded);
 }catch(e){if(e.name==='AbortError')throw new Error('参考数据下载超时，请检查网络后重试。');throw e;}finally{clearTimeout(timeout);}
}
async function cached(path){if(cache.has(path))return cache.get(path);const p=jsonFetch(path,true);cache.set(path,p);try{const value=await p;if(cache.size>16){for(const key of cache.keys()){if(!key.endsWith('index.json.gz')&&key!==path){cache.delete(key);break;}}}return value;}catch(e){cache.delete(path);throw e;}}
function transcriptNotes(r){
 if(r.manual)return r.mrna?['用户提供 mRNA；未核对版本、UTR 注释或基因组坐标。']:['尚未提供 mRNA。'];
 const m=r.mrnaInfo;if(!r.mrna){const reasons={no_exon_or_utr:'源注释缺少外显子及 UTR 信息，暂不能提供可核对的 mRNA。',not_loaded:'此参考库尚未载入 mRNA 数据。',load_error:'mRNA 读取失败，可重新查询；CDS 和上游区仍可使用。'};return [reasons[m?.reason]||'转录本结构与 CDS 无法一致核对，暂不提供该 mRNA。'];}
 const w=m?.warnings||[],notes=[];
 if(w.includes('no_5utr')&&w.includes('no_3utr'))notes.push('当前注释未给出两端 UTR，剪接序列与 CDS 相同；不代表已验证的完整 mRNA。');
 else{if(w.includes('no_5utr'))notes.push('当前模型没有注释 5′UTR，不能据此认定 5′端完整。');if(w.includes('no_3utr'))notes.push('当前模型没有注释 3′UTR，不能据此认定 3′端完整。');}
 if(w.includes('partial_cds_phase'))notes.push('CDS 首段不完整，5′非 CDS 区边界需核对。');
 if(w.includes('ambiguous_bases'))notes.push('mRNA 含不确定碱基。');return notes;
}
async function loadGene(version,id){
 const index=await cached(`${version}/index.json.gz`),tid=index[id.trim().toLowerCase()];
 if(!tid)throw new Error(`${labels[version]} 未找到编号“${id}”。请核对版本、编号或源注释是否异常；不会用其他编号代替。`);
 const settings=manifest?.transcripts?.ready&&manifest.transcripts.versions?.[version];
 const supplement=settings?(settings.available===0?Promise.resolve({status:'unavailable',reason:'no_exon_or_utr'}):cached(`${manifest.transcripts.path}/${version}/${(parseInt(shardFor(tid),16)%settings.shards).toString(16).padStart(2,'0')}.json.gz`).then(data=>{if(!data[tid])throw new Error('mRNA 索引不一致。');return data[tid];}).catch(()=>({status:'unavailable',reason:'load_error'}))):Promise.resolve(null);
 const [rows,info]=await Promise.all([cached(`${version}/${shardFor(tid)}.json.gz`),supplement]),record=rows[tid];
 if(!record)throw new Error('参考索引与序列不一致，请稍后重试。');
 try{return withTranscript({...record,version},info);}catch{return withTranscript({...record,version},{status:'unavailable',reason:'cds_sequence_mismatch'});}
}
function guess(id){if(/^Ach\d/i.test(id))return 'hy2';if(/^Actinidia\d/i.test(id))return 'hy3';if(/^(?:gene-)?(?:CEY00_)?Acc\d/i.test(id)||/Acc\d+\.\d+_mrna/.test(id))return 'red5';if(/^DTZ79_/i.test(id))return 'white';throw new Error('无法自动识别这个编号，请手动选择来源版本。');}
async function loadCatalog(){
 if(manifest)return manifest;manifest=await jsonFetch('manifest.json');const box=$('reference-catalog');box.replaceChildren();
 for(const v of manifest.versions){const p=el('div',undefined,'catalog-entry'),a=el('a',v.label);a.href=v.source;a.target='_blank';a.rel='noopener';p.append(a,document.createTextNode(manifest.ready===false?' · 自动编号查询尚未开放':` · ${v.geneCount.toLocaleString()} 个基因 / ${v.transcriptCount.toLocaleString()} 个转录本`));if(v.sequenceAudit?.supplied_mismatch)p.append(el('div',`${v.sequenceAudit.supplied_mismatch} 个模型与配套 CDS 文件有差异，对应结果会提示。`));if(v.sequenceAudit?.publishedProteinMismatch)p.append(el('div',`${v.sequenceAudit.publishedProteinMismatch} 个模型的翻译与发布蛋白有差异，对应结果会提示。`));const m=manifest.transcripts?.versions?.[v.id];if(m)p.append(el('div',`可核对 mRNA：${m.available.toLocaleString()} 条；注释不足暂不可用：${m.unavailable.toLocaleString()} 条。`));box.append(p);}
 if(manifest.ready===false){$('manual-panel').open=true;for(const id of ['gene-id','gene-version','gene-search','clone-example'])$(id).disabled=true;$('clone-intro-description').textContent='直接粘贴已有的 CDS、mRNA 和上游区序列，比较不同版本，并为选定模板生成克隆引物；自动编号查询尚未开放。';box.append(el('p','当前可使用下方的手动序列比较与克隆引物初筛。四版本参考库尚未公开接入，暂不按基因编号返回序列。'));}
 else{$('clone-intro-description').textContent='查询四个猕猴桃参考版本的同源候选，比较 CDS、mRNA 和上游区；确定模板后，一键进入载体引物设计与测序验证。';box.append(el('p',`建库日期：${manifest.built}。Red5 使用 GenBank GCA_003024255.1（PS1_1.69.0），保留 Acc 及 NCBI 编号别名。`));}
 box.append(el('p','手动导入的记录不自动确认版本、同源关系或基因组坐标。基因编号检索库使用 DIAMOND sensitive 蛋白候选，E≤1e−10、identity≥30%、双向覆盖均≥40%，保留前三名；未整合共线性。'));
 return manifest;
}
loadCatalog().then(m=>{if(!source)say('gene-status',m.ready===false?'自动编号查询尚未开放。可先在下方粘贴序列，进行比较和引物设计。':`四个版本已就绪，共 ${m.versions.reduce((n,v)=>n+v.geneCount,0).toLocaleString()} 个有效基因。`);}).catch(e=>{$('manual-panel').open=true;say('gene-status',e.message+' 也可在下方手动输入序列。',true);});
for(const id of ['gene-id','gene-version'])$(id).addEventListener(id==='gene-id'?'input':'change',()=>{clearGenes();say('gene-status','输入编号后查询四个版本。');});
$('clone-example').addEventListener('click',()=>{$('gene-id').value='Actinidia03760';$('gene-version').value='hy3';$('gene-form').dispatchEvent(new Event('submit',{cancelable:true}));});
$('gene-form').addEventListener('submit',async e=>{
 e.preventDefault();clearGenes();const rev=revision,id=$('gene-id').value.trim();if(!id)return say('gene-status','请填写一个基因或转录本编号。',true);
 $('gene-search').disabled=true;say('gene-status','正在查询参考编号与同源候选，首次需下载参考文件，请稍候…');
 try{await loadCatalog();if(manifest.ready===false)throw new Error('自动编号查询尚未开放；请先使用下方手动序列入口。');const version=$('gene-version').value==='auto'?guess(id):$('gene-version').value;const record=await loadGene(version,id);if(rev!==revision)return;source=record;const next=new Map([[version,record]]),lists=new Map();const missing=[];
 await Promise.all(VERSIONS.filter(([v])=>v!==version).map(async([v])=>{const hits=record.homologs?.[v]||[];lists.set(v,hits);if(!hits.length){missing.push(labels[v]);return;}try{const r=await loadGene(v,hits[0].id);next.set(v,{...r,evidence:hits[0]});}catch(error){missing.push(`${labels[v]}：${error.message}`);}}));
 if(rev!==revision)return;selected=next;candidateLists=lists;renderCards();renderTemplateOptions();$('gene-results').hidden=false;say('gene-status',`已读取 ${selected.size} 个版本的候选。${missing.length?'未取得候选：'+missing.join('；')+'。':''} 请核对编号与证据后比较。`);say('compare-status','选择 CDS、mRNA 或上游区，点击比较当前候选。');
 }catch(error){if(rev===revision)say('gene-status',error.message,true);}finally{$('gene-search').disabled=manifest?.ready===false;}
});
function option(value,text){const n=el('option',text);n.value=value;return n;}
function renderCards(){$('candidate-title').textContent=source.manual?'各版本输入序列':'四版本候选基因';$('candidate-note').textContent=source.manual?'手动输入不自动证明版本来源或同源关系。':'互为最佳匹配是同源候选证据，不等于已证实的直系同源。';const box=$('gene-cards');box.replaceChildren();for(const [v,label] of VERSIONS){const r=selected.get(v),card=el('article',undefined,'gene-card'+(source.version===v?' source':'')),head=el('div',undefined,'gene-card-heading');head.append(el('h3',label),el('span',r?.manual?'手动输入':source.version===v?'查询来源':'同源候选','section-number'));card.append(head);
 if(!r){card.append(el('p',source.manual?'该版本尚未提供序列。':'当前阈值下没有找到可用候选；不能据此判断该版本中不存在这个基因。'));box.append(card);continue;}
 card.append(el('strong',r.id,'gene-name'),el('p',r.manual?`用户提供序列 · CDS ${r.cds.length.toLocaleString()} bp · 上游 ${r.upstream.length.toLocaleString()} bp · 未核对基因组坐标`:`转录本：${r.transcriptLabel||r.transcript}\n${r.chrom}:${r.geneStart.toLocaleString()}–${r.geneEnd.toLocaleString()} (${r.strand}) · CDS ${r.cds.length.toLocaleString()} bp · 上游 ${r.upstream.length.toLocaleString()} bp`));
 card.append(el('p',r.mrna?`mRNA ${r.mrna.length.toLocaleString()} bp · 剪接序列 / cDNA（T）`:'mRNA 暂不可用','transcript-note'));for(const note of transcriptNotes(r))card.append(el('p',note,'transcript-note'));
 if(r.description)card.append(el('p',r.description));
 const h=r.evidence,ev=el('div',undefined,'gene-evidence');if(h){ev.append(el('span',h.reciprocal?(h.tiedBest?'双向最佳（有同分候选）':'双向最佳候选'):'单向相似候选',h.reciprocal&&!h.tiedBest?'':'caution'),el('span',`蛋白一致率 ${h.identity.toFixed(1)}%`),el('span',`覆盖 ${h.queryCoverage.toFixed(1)}% / ${h.targetCoverage.toFixed(1)}%`),el('span',`E ${h.evalue===0?'≈0':h.evalue.toExponential(1)}`));card.append(ev);}
 if(v!==source.version&&!r.manual){const list=candidateLists.get(v)||[],lab=el('label','切换候选（按蛋白匹配得分排序）'),select=el('select');select.setAttribute('aria-label',label+' 同源候选');for(const hit of list)select.append(option(hit.id,`${hit.id} · ${hit.identity.toFixed(1)}%${hit.reciprocal?' · 双向最佳':''}`));select.value=r.id;select.addEventListener('change',()=>chooseCandidate(v,select.value,list.find(x=>x.id===select.value)));lab.append(select);card.append(lab);}
 if(r.isoforms?.length>1){const lab=el('label','转录本（同源搜索使用本基因最长 CDS 蛋白）'),select=el('select');select.setAttribute('aria-label',label+' 转录本');for(const tid of r.isoforms)select.append(option(tid,tid.split('|').at(-1).replace(/_mrna$/,'')+(tid===r.representative?' · 最长 CDS':'')));select.value=r.transcript;select.addEventListener('change',()=>chooseCandidate(v,select.value,r.evidence));lab.append(select);card.append(lab);}
 for(const w of r.warnings||[])card.append(el('p',w,'gene-warning'));box.append(card);
 }}
async function chooseCandidate(v,id,evidence){invalidateComparison();const rev=revision;say('gene-status',`正在读取 ${id}…`);try{const r=await loadGene(v,id);if(rev!==revision)return;selected.set(v,{...r,evidence});renderCards();renderTemplateOptions();say('gene-status','候选已更新，请重新比较并设计引物。');}catch(e){if(rev===revision){say('gene-status',e.message,true);renderCards();}}}
function kind(){return $('clone-kind').value;}
function upstreamLength(){return Number($('promoter-length').value);}
function rows(){return VERSIONS.map(([v])=>selected.get(v)).filter(Boolean);}
function usableRows(){return rows().filter(r=>kind()==='cds'?r.cds:kind()==='mrna'?r.mrna:r.upstream);}
function currentTemplate(){const r=selected.get($('clone-template').value);if(!r)throw new Error('请先查询并选择一个模板。');return {record:r,...cloningTemplate(r,kind(),upstreamLength(),$('remove-stop').checked)};}
function renderTemplateOptions(){const old=$('clone-template').value,box=$('clone-template'),list=usableRows();box.replaceChildren();for(const r of list)box.append(option(r.version,`${labels[r.version]} · ${r.id} · ${r.transcriptLabel||r.transcript}`));box.value=list.some(r=>r.version===old)?old:list.some(r=>r.version==='white')?'white':list[0]?.version||'';updateTemplate();}
function updateTemplate(){
 try{
  const t=currentTemplate(),r=t.record;
  const material=kind()==='promoter'?'请使用基因组 DNA 克隆上游区':kind()==='mrna'?'mRNA 对应 cDNA 克隆模板':'请使用 cDNA 克隆 CDS';
  $('template-summary').textContent=`${labels[r.version]} / ${r.id} · ${t.sequence.length.toLocaleString()} bp · ${material}`+(t.removedStop?' · 已去除末端终止密码子':'')+(t.truncated?` · 边界限制：请求 ${t.requested} bp，实际仅 ${t.sequence.length} bp`:'')+(kind()==='promoter'?(t.start===null?' · 用户提供上游区，坐标未核对':` · ${r.chrom}:${t.start}–${t.end} (${r.strand})`):'')+(r.warnings?.length?' · 注意：'+r.warnings.join('；'):'')+(kind()==='mrna'?' · '+transcriptNotes(r).join('；'):'')+(r.evidence&&!r.evidence.reciprocal?' · 当前为单向蛋白相似候选，请结合基因结构或共线性确认克隆对象。':'');
  renderSequence(t);$('design-primers').disabled=false;
 }catch(e){$('template-summary').textContent=kind()==='mrna'?'当前候选没有可用的 mRNA。可切换候选，或在手动输入区提供已核对的转录本序列。':e.message;$('template-sequence').replaceChildren();$('sequence-legend').replaceChildren();$('sequence-view-note').textContent='';$('design-primers').disabled=true;}
}
function templateFasta(t){return `>${t.record.version}|${t.record.id}|${t.record.transcript}|${t.kind}|${t.sequence.length}bp${t.kind==='promoter'&&t.start!==null?`|${t.record.chrom}:${t.start}-${t.end}(${t.record.strand})|upstream_only`:t.record.manual?'|user_provided':''}${t.kind==='mrna'?`|spliced_cDNA|CDS=${t.cdsStart??'unknown'}..${t.cdsEnd??'unknown'}`:''}${t.removedStop?'|stop_removed':''}\n${t.sequence.match(/.{1,80}/g)?.join('\n')||''}`;}
function renderSequence(t){
 const r=t.record,view=sequenceDisplay(r,t.kind,upstreamLength(),$('remove-stop').checked),pre=$('template-sequence'),legend=$('sequence-legend');pre.replaceChildren();legend.replaceChildren();
 let header=templateFasta(t).split('\n')[0],note='';
 if(t.kind==='promoter'){
  header=`>${r.version}|${r.id}|${r.transcript}|UPSTREAM_PLUS_SPLICED_CDS_VIEW|upstream=1..${view.upstreamLength}|CDS=${view.cdsStart??'absent'}..${view.cdsEnd??'absent'}|not_genomic_contig`;
  note=`显示上游 ${view.upstreamLength.toLocaleString()} bp`+(r.cds?`＋拼接 CDS ${r.cds.length.toLocaleString()} bp（红色）。两段拼接供对照，不代表连续基因组片段。`:'；此记录未提供 CDS。')+`当前比较、引物和模板 FASTA 仅使用上游 ${t.sequence.length.toLocaleString()} bp。`;
  legend.append(el('span',`上游区：1–${view.upstreamLength} bp`));
 }else if(t.kind==='mrna'){
  note='按基因转录方向显示剪接 mRNA / cDNA，使用 T；不含内含子。不自动添加未注释的 UTR 或 poly(A) 尾。';
  if(view.cdsStart>1)legend.append(el('span',`5′非 CDS 区：1–${view.cdsStart-1} bp`));
  if(!view.cdsStart)note+=' 当前未核定 CDS 在这条 mRNA 内的唯一位置。';
 }else note='当前显示与克隆模板一致的 CDS。';
 if(view.cdsStart&&view.cdsEnd){legend.append(el('span',`CDS（红色）：${view.cdsStart}–${view.cdsEnd} bp`,'cds-red'));}
 if(t.kind==='mrna'&&view.cdsEnd&&view.cdsEnd<view.sequence.length)legend.append(el('span',`3′非 CDS 区：${view.cdsEnd+1}–${view.sequence.length} bp`));
 $('sequence-view-note').textContent=note;pre.setAttribute('aria-label',note+(view.cdsStart?` CDS 从第 ${view.cdsStart} 至 ${view.cdsEnd} 位。`:''));
 pre.append(document.createTextNode(header+'\n'));
 const start=view.cdsStart?view.cdsStart-1:-1,end=view.cdsEnd??-1;
 for(let i=0;i<view.sequence.length;i+=80){
  const finish=Math.min(i+80,view.sequence.length),cuts=[i,finish];if(start>i&&start<finish)cuts.push(start);if(end>i&&end<finish)cuts.push(end);cuts.sort((a,b)=>a-b);
  for(let j=0;j<cuts.length-1;j++){const a=cuts[j],b=cuts[j+1],text=view.sequence.slice(a,b);pre.append(a>=start&&a<end?el('span',text,'cds-red'):document.createTextNode(text));}
  if(finish<view.sequence.length)pre.append(document.createTextNode('\n'));
 }
}
$('clone-kind').addEventListener('change',()=>{
 invalidateComparison();const promoter=kind()==='promoter',mrna=kind()==='mrna';$('promoter-setting').hidden=!promoter;document.querySelector('.stop-setting').hidden=kind()!=='cds';if(kind()!=='cds')$('remove-stop').checked=false;
 $('sequence-definition').textContent=promoter?'比较注释 CDS 起点之前的上游区，长度上限 2000 bp。下方查看区追加红色 CDS 作为边界对照，不计入启动子相似率或引物产物。':mrna?'mRNA 按外显子或 CDS＋UTR 注释拼接，含已注释的非编码区，去除内含子；采用 T 表示 cDNA。只有可核对的转录本参加比较，缺少 UTR 不代表两端完整。':'CDS 从同版本基因组按注释拼接，并统一到转录方向。全长一致率将缺口和不确定碱基计入分母。';
 if(source?.manual)$('sequence-definition').textContent=promoter?'使用您提供的上游区，按 5′→3′ 比较；缩短时保留末端目标长度。查看区追加您提供的红色 CDS，不改变上游比较/引物范围。':mrna?'使用您提供的 mRNA / cDNA，按输入的 5′→3′ 方向比较。未自动核对 UTR、转录本结构或同源关系。':'使用您提供的 CDS，未自动验证版本、注释或同源关系。';
 renderTemplateOptions();if(promoter||mrna)$('template-detail').open=true;
});
$('promoter-length').addEventListener('input',()=>{invalidateComparison();updateTemplate();});
$('clone-template').addEventListener('change',()=>{revision++;cancelJobs();clearPrimer();updateTemplate();});
$('remove-stop').addEventListener('change',()=>{revision++;cancelJobs();clearPrimer();updateTemplate();});
function task(type,data){
 if(typeof Worker==='undefined')return new Promise((resolve,reject)=>setTimeout(()=>{try{if(type==='compare')resolve(compareDNA(data.a,data.b));else{const r=designPrimers(data.sequence,data.options);for(const p of r.pairs)p.binding=data.references.map(x=>({name:x.name,forward:bindingCheck(p.forward.sequence,x.sequence),reverse:bindingCheck(p.reverse.sequence,x.sequence,true)}));resolve(r);}}catch(e){reject(e);}},0));
 if(!worker){worker=new Worker(new URL('./clone-worker.mjs?v=20261019-clone4',import.meta.url),{type:'module'});worker.onmessage=({data})=>{const j=jobs.get(data.id);if(!j)return;clearTimeout(j.timer);jobs.delete(data.id);data.error?j.reject(new Error(data.error)):j.resolve(data.result);};worker.onerror=()=>{for(const j of jobs.values()){clearTimeout(j.timer);j.reject(new Error('后台计算无法启动，请刷新页面重试。'));}jobs.clear();worker.terminate();worker=null;};}
 return new Promise((resolve,reject)=>{const id=++jobId;const timer=setTimeout(()=>{jobs.delete(id);reject(new Error('计算超时，请缩短片段或下载后离线处理。'));},120000);jobs.set(id,{resolve,reject,timer});worker.postMessage({id,type,...data});});
}
$('compare-genes').addEventListener('click',async()=>{
 invalidateComparison();const rev=revision;$('compare-genes').disabled=true;const list=usableRows(),out=[];
 try{if(list.length<2)throw new Error('至少需要两个版本的有效候选。');const seqs=list.map(r=>cloningTemplate(r,kind(),upstreamLength()).sequence);let done=0,total=list.length*(list.length-1)/2;
 for(let i=0;i<list.length;i++)for(let j=i+1;j<list.length;j++){say('compare-status',`正在比较 ${labels[list[i].version]} × ${labels[list[j].version]}（${++done}/${total}）…`);try{const result=await task('compare',{a:seqs[i],b:seqs[j]});if(rev!==revision)return;out.push({i,j,a:list[i],b:list[j],result});}catch(e){if(rev!==revision)return;out.push({i,j,a:list[i],b:list[j],error:e.message});}}
 if(rev!==revision)return;comparisons=out;renderMatrix(list);say('compare-status',`已完成 ${out.filter(x=>x.result).length} / ${total} 组全长比较。点击分数查看逐位差异。${out.some(x=>x.error)?' 部分序列超出计算范围，未生成分数。':''}`);
 }catch(e){if(rev===revision)say('compare-status',e.message,true);}finally{$('compare-genes').disabled=false;}
});
function renderMatrix(list){const table=el('table'),caption=el('caption',`${kind()==='cds'?'CDS':kind()==='mrna'?'mRNA / cDNA':'CDS 上游 '+upstreamLength()+' bp'} · 全长核苷酸一致率（缺口计入分母）`),head=el('thead'),tr=el('tr');tr.append(el('th','参考版本'));for(const r of list)tr.append(el('th',labels[r.version]));head.append(tr);table.append(caption,head);const body=el('tbody');for(let i=0;i<list.length;i++){const row=el('tr'),th=el('th',labels[list[i].version]);th.scope='row';row.append(th);for(let j=0;j<list.length;j++){const td=el('td');if(i===j)td.textContent='—';else{const p=comparisons.find(p=>(p.i===i&&p.j===j)||(p.i===j&&p.j===i));if(p?.result){const r=p.result,b=el('button',undefined,'matrix-cell'+(r.identity<70?' low':r.identity<90?' medium':''));b.type='button';b.append(el('strong',r.identity.toFixed(2)+'%'),el('small','无缺口覆盖 '+Math.min(r.coverageA,r.coverageB).toFixed(1)+'%'));b.title='点击查看逐位比对';b.addEventListener('click',()=>showPair(p));td.append(b);}else{td.textContent='未计算';td.title=p?.error||'';}}row.append(td);}body.append(row);}table.append(body);$('comparison-matrix').replaceChildren(table);}
function showPair(p){activePair=p;const r=p.result;$('pair-detail').hidden=false;$('pair-title').textContent=`${p.a.id} × ${p.b.id}`;$('pair-stats').textContent=`全长一致率 ${r.identity.toFixed(2)}% · 相同 ${r.matches} / ${r.columns} 列 · 替换 ${r.mismatches} · 缺口 ${r.gaps} · 模糊 ${r.ambiguous} · 无缺口覆盖 ${r.coverageA.toFixed(1)}% / ${r.coverageB.toFixed(1)}%`;
 const lines=[];let ca=0,cb=0;for(let i=0;i<r.columns;i+=60){const a=r.alignedA.slice(i,i+60),b=r.alignedB.slice(i,i+60),sa=ca+1,sb=cb+1;ca+=a.replace(/-/g,'').length;cb+=b.replace(/-/g,'').length;lines.push(`${labels[p.a.version].padEnd(13)} ${String(sa).padStart(5)} ${a} ${ca}`,`                    ${[...a].map((x,j)=>x===b[j]&&'ACGT'.includes(x)?'|':x==='-'||b[j]==='-'?' ':'·').join('')}`,`${labels[p.b.version].padEnd(13)} ${String(sb).padStart(5)} ${b} ${cb}`,'');} $('pair-alignment').textContent=lines.join('\n');}
function save(name,text,type='text/plain;charset=utf-8'){const url=URL.createObjectURL(new Blob([text],{type})),a=el('a');a.href=url;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),2000);}
function checkedAction(id,fn){$(id).addEventListener('click',()=>{try{fn();}catch(e){say('compare-status',e.message,true);}});}
checkedAction('clone-download-fasta',()=>save('KiwiGene_four_versions_'+kind()+'.fasta',fastaText(usableRows(),kind(),upstreamLength())));
checkedAction('template-download',()=>{const t=currentTemplate();save(`KiwiGene_${t.record.id}_${t.kind}.fasta`,templateFasta(t));});
checkedAction('annotated-download',()=>{
 const t=currentTemplate(),article=el('article');article.append(el('h1',`${labels[t.record.version]} · ${t.record.id}`),el('p',$('template-summary').textContent),el('p',$('sequence-view-note').textContent),$('sequence-legend').cloneNode(true),$('template-sequence').cloneNode(true));
 const html='<!DOCTYPE html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>猕猴桃基因工坊 · 标注序列</title><style>body{max-width:1000px;margin:30px auto;padding:0 20px;font:16px/1.8 system-ui;color:#253b31}h1{font-size:24px}pre{font:13px/1.8 monospace;white-space:pre-wrap;overflow-wrap:anywhere;padding:20px;background:#f7f9f5;border:1px solid #dce5d9}.cds-red{color:#b42318;font-weight:bold}.sequence-legend{display:flex;gap:12px;flex-wrap:wrap}.sequence-legend span{padding:6px;border:1px solid #ddd}@media print{.cds-red{color:#b42318!important;print-color-adjust:exact;-webkit-print-color-adjust:exact}}</style><body>'+article.outerHTML+'</body></html>';
 save(`KiwiGene_${t.record.id}_${t.kind}_annotated.html`,html,'text/html;charset=utf-8');
});
checkedAction('pair-download',()=>{if(!activePair)throw new Error('请先选择一组比较。');const p=activePair;save(`KiwiGene_${p.a.id}_${p.b.id}_aligned.fasta`,`>${p.a.version}|${p.a.id}\n${p.result.alignedA}\n>${p.b.version}|${p.b.id}\n${p.result.alignedB}\n`);});
checkedAction('send-to-vector',()=>{const record=selected.get($('clone-template').value);window.dispatchEvent(new CustomEvent('seqlab-vector-input',{detail:vectorInputFromRecord(record,kind(),upstreamLength())}));});
checkedAction('send-to-alignment',()=>{const t=currentTemplate();window.dispatchEvent(new CustomEvent('seqlab-reference',{detail:{name:`${t.record.version}|${t.record.id}|${t.kind}`,sequence:t.sequence}}));window.dispatchEvent(new CustomEvent('seqlab-show-workspace',{detail:{workspace:'sequence'}}));$('sequence-input').scrollIntoView({behavior:'smooth',block:'center'});});
for(const id of ['primer-tm','primer-min','primer-max','primer-salt','primer-mg','primer-dntp','primer-conc','forward-tail','reverse-tail'])$(id).addEventListener('input',()=>{revision++;cancelJobs();clearPrimer();});
$('design-primers').addEventListener('click',async()=>{
 clearPrimer();const rev=revision;$('design-primers').disabled=true;say('primer-status','正在筛选覆盖目标片段两端的引物，并检查当前候选序列中的匹配…');
 try{const t=currentTemplate(),options={target:Number($('primer-tm').value),minLen:Number($('primer-min').value),maxLen:Number($('primer-max').value),salt:Number($('primer-salt').value),mg:Number($('primer-mg').value),dntp:Number($('primer-dntp').value),concentration:Number($('primer-conc').value),forwardTail:$('forward-tail').value,reverseTail:$('reverse-tail').value};
 const references=usableRows().map(r=>({name:`${labels[r.version]} · ${r.id}`,sequence:cloningTemplate(r,kind(),upstreamLength()).sequence}));const result=await task('design',{sequence:t.sequence,options,references});if(rev!==revision)return;primerResult={...result,template:t};renderPrimers();if(!result.pairs.length)return say('primer-status',result.reason,true);say('primer-status',`已生成 ${result.pairs.length} 组初筛候选。尚未进行全基因组特异性验证。${result.templateAmbiguous?' 模板含 '+result.templateAmbiguous+' 个不确定碱基，需核对实测序列。':''}`);
 }catch(e){if(rev===revision)say('primer-status',e.message,true);}finally{$('design-primers').disabled=false;}
});
function bindingText(x){if(!x)return '未找到可比较区间';return `${x.start}–${x.end}，错配 ${x.mismatches}，3′末 5 nt 错配 ${x.endMismatches}${x.matches>1?'，'+x.matches+' 处同分位置':''}`;}
function renderPrimers(){const result=primerResult,t=result.template,box=$('primer-results');box.replaceChildren();const lines=['猕猴桃基因工坊 · 克隆引物初筛报告',`${labels[t.record.version]} | ${t.record.id} | ${t.record.transcript}`,`对象：${t.kind}；长度：${t.sequence.length} bp；方向：5′→3′`,kind()==='promoter'?(t.start===null?'用户提供上游区，未核对基因组坐标。':`坐标：${t.record.chrom}:${t.start}-${t.end} (${t.record.strand})；注释 CDS 上游，不是已验证 TSS 启动子。`):`${t.kind==='mrna'?'mRNA':'CDS'} 模板：cDNA；`+(t.removedStop?'已去终止密码子。':'保留原始末端。'),`参数：${JSON.stringify(result.parameters)}`,'方法：最近邻 Tm；连续互补 / 简化发夹初筛。不是 Primer3 输出。','未进行完整二级结构热力学或全基因组特异性验证；以下匹配仅在选中的候选序列中搜索。',...(t.record.warnings||[]).map(w=>'模板警示：'+w),...(t.kind==='mrna'?transcriptNotes(t.record):[]).map(w=>'mRNA 提示：'+w),''];
 result.pairs.forEach((p,i)=>{const card=el('article',undefined,'primer-pair');card.append(el('h3',`候选 ${i+1} · 目标片段 ${p.productSize} bp${p.productWithTails!==p.productSize?' · 含附加序列 '+p.productWithTails+' bp':''}`));lines.push(`候选 ${i+1} | 目标 ${p.productSize} bp | 含附加序列 ${p.productWithTails} bp`);
 for(const [dir,pr] of [['F 正向',p.forward],['R 反向',p.reverse]]){const row=el('div',undefined,'primer-row'),detail=el('div');detail.append(el('code',pr.fullSequence));const stats=`退火区 ${pr.length} nt · Tm ${pr.tm.toFixed(2)} °C · GC ${pr.gc.toFixed(1)}% · 3′自互补 ${pr.self.end} bp · 连续发夹茎 ${pr.hairpin} bp`;detail.append(el('small',stats));if(pr.tail)detail.append(el('small',`5′附加：${pr.tail}；退火区：${pr.sequence}`));row.append(el('span',dir+" 5′→3′"),detail);card.append(row);lines.push(dir+': '+pr.fullSequence,stats,'退火区：'+pr.sequence);}
 card.append(el('p',`ΔTm ${p.tmDifference.toFixed(2)} °C · 两引物 3′ 连续互补最长 ${p.heterodimer.end} bp · 简化结构筛查通过`,'field-hint'));
 const bindings=el('div',undefined,'binding-report');bindings.append(el('strong','在当前候选序列中匹配（仅初筛）'));for(const b of p.binding){const text=`${b.name}：F ${bindingText(b.forward)}；R ${bindingText(b.reverse)}`;bindings.append(el('p',text));lines.push(text);}card.append(bindings);box.append(card);lines.push('');});
 lines.push('模板 FASTA：',templateFasta(t));primerReport=lines.join('\n');$('primer-downloads').hidden=!result.pairs.length;
}
$('primer-copy').addEventListener('click',async()=>{try{if(!primerReport)return;await navigator.clipboard.writeText(primerReport);say('primer-status','引物及报告已复制。');}catch{say('primer-status','浏览器未允许复制，请点击下载引物报告。',true);}});
$('primer-download').addEventListener('click',()=>{if(primerReport)save('KiwiGene_'+primerResult.template.record.id+'_primers.txt',primerReport);});
$('primer-tsv').addEventListener('click',()=>{if(!primerResult?.pairs.length)return;const rows=[['名称','序列_5to3','退火区','Tm_C','GC_percent','目标产物_bp','参考版本','基因编号','验证状态']];primerResult.pairs.forEach((p,i)=>{for(const [d,x] of [['F',p.forward],['R',p.reverse]])rows.push([`${primerResult.template.record.id}_${i+1}_${d}`,x.fullSequence,x.sequence,x.tm.toFixed(2),x.gc.toFixed(1),p.productSize,primerResult.template.record.version,primerResult.template.record.id,'未做全基因组特异性验证']);});save('KiwiGene_primer_order.tsv','\uFEFF'+rows.map(r=>r.join('\t')).join('\n'),'text/tab-separated-values;charset=utf-8');});

let simulatedInput=false;
for(const [v,name] of VERSIONS){
 const box=el('div',undefined,'manual-version');box.append(el('h3',name));
 for(const [key,title,tag] of [['id','基因名称／编号（可选）','input'],['cds','CDS 序列','textarea'],['mrna','mRNA / cDNA 序列（含已知 UTR，使用 T）','textarea'],['upstream','上游区序列（基因转录方向，≤2000 bp）','textarea']]){const label=el('label',title),input=el(tag);input.id=`manual-${v}-${key}`;input.autocomplete='off';input.spellcheck=false;input.setAttribute('aria-label',name+' '+title);if(tag==='textarea')input.placeholder='>序列名称\nACGT…';input.addEventListener('input',()=>{simulatedInput=false;$('manual-demo-note').hidden=true;if(source?.manual){clearGenes();say('manual-status','序列已修改，请重新点击“使用这些序列”。');}});label.append(input);box.append(label);}
 $('manual-fields').append(box);
}
function manualSequence(v,key){const input=$(`manual-${v}-${key}`).value.trim();if(!input)return '';const title={cds:'CDS',mrna:'mRNA',upstream:'上游区'}[key],parsed=parseSequences(input);if(parsed.length!==1)throw new Error(`${labels[v]} 的 ${title} 每格只接收一条序列。`);const seq=parsed[0].sequence,max=key==='upstream'?2000:12000;if(seq.length<60||seq.length>max)throw new Error(`${labels[v]} 的 ${title} 需要 60–${max===2000?'2000':'12,000'} bp，当前 ${seq.length} bp。`);return seq;}
$('manual-use').addEventListener('click',()=>{
 clearGenes();try{const records=[];for(const [v] of VERSIONS){const cds=manualSequence(v,'cds'),mrna=manualSequence(v,'mrna'),upstream=manualSequence(v,'upstream');if(!cds&&!mrna&&!upstream)continue;const id=$(`manual-${v}-id`).value.trim().replace(/[\s>|]/g,'_').slice(0,80)||'USER_'+v;const warnings=['用户提供序列；尚未核对参考版本、坐标或同源关系。'];if(simulatedInput)warnings.unshift('模拟序列，不代表真实基因。');if(cds&&(!cds.startsWith('ATG')||cds.length%3||!/T(?:AA|AG|GA)$/.test(cds)))warnings.push('CDS 起止密码子或长度不符合完整开放阅读框，请核对。');const pos=cds&&mrna?mrna.indexOf(cds):-1;const mrnaInfo=pos>=0&&pos===mrna.lastIndexOf(cds)?{status:'available',cdsStart:pos+1,cdsEnd:pos+cds.length}:{status:mrna?'available':'unavailable',cdsStart:null,cdsEnd:null};records.push({id,transcript:id,version:v,manual:true,cds,mrna,mrnaInfo,upstream,strand:'+',chrom:null,geneStart:null,geneEnd:null,upstreamStart:null,upstreamEnd:null,warnings,isoforms:[id],representative:id});}
 if(!records.length)throw new Error('请至少粘贴一个版本的一条 CDS、mRNA 或上游区序列。');selected=new Map(records.map(r=>[r.version,r]));source=records[0];candidateLists=new Map();$('clone-kind').value=records.some(r=>r.cds)?'cds':records.some(r=>r.mrna)?'mrna':'promoter';renderCards();renderTemplateOptions();$('clone-kind').dispatchEvent(new Event('change'));$('gene-results').hidden=false;say('gene-status',`正在使用 ${records.length} 个手动输入的版本；没有运行自动同源检索。`);say('manual-status',`已载入 ${records.length} 个版本。至少两条可比较；单条即可设计引物。`);say('compare-status','选择比较对象后开始比较，或直接生成当前模板的引物。');$('manual-demo-note').hidden=!simulatedInput;
 }catch(e){say('manual-status',e.message,true);}
});
$('manual-example').addEventListener('click',()=>{$('manual-panel').open=true;
 const codons=['GCC','GCT','CAG','ATC','GAC','GGT','CTG','AAC'];let seed=12073;const bases='ACGT';let up='';for(let i=0;i<2000;i++){seed=(Math.imul(seed,1664525)+1013904223)>>>0;up+=bases[seed>>>30];}
 VERSIONS.forEach(([v],i)=>{const middle=Array.from({length:189},(_,k)=>codons[k%codons.length]);middle[20+i*10]=codons[(i+3)%codons.length];$(`manual-${v}-id`).value='SIMULATED_'+v.toUpperCase();$(`manual-${v}-cds`).value='ATGACCGTACAGGATCTGCAG'+middle.join('')+['TGC','ATC','AGG','TAC','CGA','GCT','TAA'].join('');$(`manual-${v}-mrna`).value='CTAGCTAGGATCAGTCTAGG'+$(`manual-${v}-cds`).value+'TTCAGTACGATCTGACTGAT';const at=500+i*73;$(`manual-${v}-upstream`).value=up.slice(0,at)+bases[(bases.indexOf(up[at])+1)%4]+up.slice(at+1);});simulatedInput=true;$('manual-use').click();
});

