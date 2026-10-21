// Synthetic, deterministic demonstration only; not a reference kiwifruit gene.
export function vectorDemo(){
 let seed=1;
 const randomIndex=n=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return Math.floor(seed/4294967296*n);};
 const codons=['GCT','GAA','CAG','ATC','GAC','GGT','CTG','AAC','TTC','AAG','ACT','TGT','CAT','ATT','CCA','AGT','GTT','TAC','GGA','CGT'];
 const cds='ATG'+Array.from({length:110},()=>codons[randomIndex(codons.length)]).join('')+'TAA';
 const promoter=Array.from({length:540},()=>'ACGT'[randomIndex(4)]).join('');
 return {geneName:'AcDemo_simulated',cds,promoter};
}
