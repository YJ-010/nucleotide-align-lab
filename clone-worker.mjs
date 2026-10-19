import {compareDNA,designPrimers,bindingCheck} from './clone-core.mjs?v=20261019-clone3';
self.onmessage=({data})=>{const {id,type}=data;try{
 let result;
 if(type==='compare')result=compareDNA(data.a,data.b);
 else if(type==='design'){result=designPrimers(data.sequence,data.options);result.pairs.forEach(pair=>{pair.binding=(data.references||[]).map(r=>({name:r.name,forward:bindingCheck(pair.forward.sequence,r.sequence),reverse:bindingCheck(pair.reverse.sequence,r.sequence,true)}));});}
 else throw new Error('未知计算任务');
 self.postMessage({id,result});
}catch(e){self.postMessage({id,error:e.message});}};
