import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';

test('static DOM IDs match all literal controller references, including catalog description',()=>{
 const read=path=>readFileSync(new URL('../'+path,import.meta.url),'utf8');
 const html=read('index.html');
 const ids=[...html.matchAll(/\bid="([^"]+)"/g)].map(x=>x[1]);
 assert.equal(new Set(ids).size,ids.length,'HTML IDs must be unique');
 for(const name of ['app.mjs','clone-app.mjs','vector-app.mjs','workspace-nav.mjs']){
  const js=read(name);
  for(const [,id] of js.matchAll(/\$\(['"]([^'"]+)['"]\)/g))assert.ok(ids.includes(id),name+' needs #'+id);
 }
 assert.ok(ids.includes('clone-intro-description'));
 assert.ok(!read('clone-app.mjs').includes('.clone-intro > p'));
});
