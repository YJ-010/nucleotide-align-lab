import {compareDNA} from './clone-core.mjs?v=20261026-boundary1';
import {designCloningSchemes} from './cloning-design.mjs?v=20261026-boundary1';
self.onmessage=({data})=>{const {id,type}=data;try{
 let result;
 if(type==='compare')result=compareDNA(data.a,data.b);
 else if(type==='design')result=designCloningSchemes(data.record,data.kind,data.length,data.options,data.override);
 else throw new Error('未知计算任务');
 self.postMessage({id,result});
}catch(e){self.postMessage({id,error:e.message});}};
