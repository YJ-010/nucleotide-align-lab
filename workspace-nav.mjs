const WORKSPACES=['clone','vector','sequence'];
const hashes={clone:'#clone',vector:'#vector',sequence:'#alignment'};
const byHash={'#clone':'clone','#vector':'vector','#construct':'vector','#alignment':'sequence','#sequence':'sequence'};
const $=id=>document.getElementById(id);

export function showWorkspace(which,{updateHash=true,focus=false}={}){
 const active=WORKSPACES.includes(which)?which:'clone';
 for(const name of WORKSPACES){
  const panel=$(name+'-workspace'),tab=$(name+'-tab'),on=name===active;
  if(panel)panel.hidden=!on;
  if(tab){
   tab.classList.toggle('active',on);
   tab.setAttribute('aria-selected',String(on));
   tab.tabIndex=on?0:-1;
  }
 }
 if(updateHash&&window.location.hash!==hashes[active])history.replaceState(null,'',hashes[active]);
 if(focus)$(active+'-tab')?.focus();
 document.title='猕猴桃基因工坊 · '+(active==='clone'?'基因克隆':active==='vector'?'载体构建与引物设计':'测序比对验证');
 return active;
}

for(const name of WORKSPACES){
 $(name+'-tab')?.addEventListener('click',()=>showWorkspace(name));
 $(name+'-tab')?.addEventListener('keydown',event=>{
  if(!['ArrowLeft','ArrowRight','Home','End'].includes(event.key))return;
  event.preventDefault();
  const current=WORKSPACES.indexOf(name);
  const next=event.key==='Home'?0:event.key==='End'?WORKSPACES.length-1:event.key==='ArrowLeft'?(current-1+WORKSPACES.length)%WORKSPACES.length:(current+1)%WORKSPACES.length;
  showWorkspace(WORKSPACES[next],{focus:true});
 });
}

window.addEventListener('hashchange',()=>showWorkspace(byHash[window.location.hash]||'clone',{updateHash:false}));
window.addEventListener('seqlab-show-workspace',event=>showWorkspace(event.detail?.workspace||'clone'));
document.querySelectorAll('[data-open-workspace]').forEach(button=>button.addEventListener('click',()=>{showWorkspace(button.dataset.openWorkspace);document.querySelector('.workspace-nav').scrollIntoView({behavior:'smooth',block:'start'});}));
document.querySelector('.help-link')?.addEventListener('click',event=>{
 event.preventDefault();
 const active=WORKSPACES.find(name=>!$(name+'-workspace')?.hidden)||'clone';
 $({clone:'clone-help',vector:'vector-help',sequence:'help'}[active])?.scrollIntoView({behavior:'smooth',block:'start'});
});
showWorkspace(byHash[window.location.hash]||'clone',{updateHash:!byHash[window.location.hash]});

