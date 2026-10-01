/* Standalone result rendering preserves search and interpretation evidence. */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const context = vm.createContext({
  I18N: {ja: {}, en: {}}, document: {addEventListener() {}},
  renderCheck(data) { return Object.values(data.primers).join(' '); },
  renderDesign() { return ''; }, renderSequence() { return ''; }, renderAssay() { return ''; },
  runMode() {},
  t: key => key, esc: value => String(value ?? '').replaceAll('&', '&amp;').replaceAll('<', '&lt;'),
  dl: (type, name, data) => `<download type="${type}">${String(data).replaceAll('<', '&lt;')}</download>`,
  primerRows: primers => Object.values(primers).join(' ')
});
for (const file of ['locus.js', 'sequence-view.js', 'sequence-tools.js']) {
  vm.runInContext(fs.readFileSync('primerblast_oss/webapp/static/' + file, 'utf8'), context);
}
const run = code => vm.runInContext(code, context);
const blast = run(`renderBlast({task:'blastn',queries:[{id:'query1',name:'<query>'}],results:[{db:'ref',tsv:'query',at_target_limit:['query1'],hits:[{qseqid:'query1',sseqid:'chr1',sstart:80,send:20,strand:'-',pident:100,length:61,evalue:1e-8,bitscore:50,qstart:1,qend:61,qseq:'ACGT',sseq:'ACGT'}]}]})`);
assert.ok(blast.includes('80–20 (-)'));
assert.ok(blast.includes('tools.limited'));
assert.ok(blast.includes('&lt;query>'));
const product = run(String.raw`renderCheck({input_orientation:'auto',primers:{F:'AAAA',R:'AAAA'},results:[{oligos:{F:'ACGT',R:'CGTA'},reverse_complemented_inputs:['R'],search_completeness:'incomplete',fasta:'>chr1\nACGT',products:[{subject:'chr1',start:1,end:4,size:4,orientation:'F/R',fwd_primer:'F',rev_primer:'R',sequence:'ACGT',fasta:'>chr1\nACGT',annotations:{status:'unavailable',genes:[]}}]}]})`);
assert.ok(product.includes('ACGT CGTA'));
assert.ok(!product.replace(/<download\b[\s\S]*?<\/download>/g, '').includes('AAAA'));
assert.ok(product.includes('incomplete'));
assert.ok(product.includes('tools.hypothesis: R'));
assert.ok(product.includes('tools.annotationMissing'));
assert.ok(!product.includes('seq.noGenes'));
const primer3 = run(String.raw`renderPrimer3({templates:[{template_id:'template',template_sequence:'ACGTACGT',pairs:[{index:0,forward:'AC',reverse:'AC',left_start:0,left_len:2,right_start:7,right_len:2,product_size:8,tm_f:60,tm_r:60,gc_f:50,gc_r:50,penalty:1,fasta:'>pair\nACGTACGT'}]}]})`);
assert.ok(primer3.includes('tools.unscreened'));
assert.ok(primer3.includes('data-primer3-check="0"'));
assert.ok(primer3.includes('type="fasta"'));
console.log('Standalone BLAST/Primer3, orientation hypotheses and FASTA evidence passed.');
