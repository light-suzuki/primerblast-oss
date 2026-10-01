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
const swapped = run(`productInputEvidence({start:1,end:100,fwd_primer:'R',rev_primer:'F',input_evidence:{sequence_status:'as_supplied',changed_primers:[],label_orientation:'swapped'},primer_bindings:[{primer:'R',input_sequence:'TTGC',oligo:'TTGC',reference_strand:'+',end5:1,end3:4,mismatches:0},{primer:'F',input_sequence:'ACGA',oligo:'ACGA',reference_strand:'-',end5:100,end3:97,mismatches:0}]})`);
assert.ok(swapped.includes('tools.literalProduct') && swapped.includes('tools.swapped'));
assert.ok(swapped.includes('R ＋ →') && swapped.includes('← − F'));
assert.ok(!swapped.includes('tools.usedSequence'));
const corrected = run(`productInputEvidence({start:1,end:100,fwd_primer:'F',rev_primer:'R',input_evidence:{sequence_status:'reverse_complement_candidate',changed_primers:['R'],label_orientation:'as_labeled',original_search_complete:false,as_supplied_locus_observed:false},primer_bindings:[{primer:'R',input_sequence:'TTGC',oligo:'GCAA',reference_strand:'-',end5:100,end3:97,mismatches:0}]})`);
assert.ok(corrected.includes('tools.changedProduct') && corrected.includes('tools.unknownOriginalLocus'));
assert.ok(corrected.includes('TTGC') && corrected.includes('GCAA') && corrected.includes('tools.usedSequence'));
assert.ok(!corrected.includes('tools.noOriginalLocus'));
const summary = run(`inputAssessments({input_assessments:[{db:'ref',as_supplied_products:0,distinct_primer_products:0,same_primer_products:0,reverse_complement_candidates:1,original_search_complete:false}]})`);
assert.ok(summary.includes('tools.unknownLiteral') && summary.includes('tools.errorUnknown'));
assert.ok(!summary.includes('tools.noLiteral'));
const sites = run(`bindingSiteEvidence({products:[],binding_site_counts:{F:{'+':0,'-':1}},binding_sites:[{primer:'F',subject:'<chr1>',reference_strand:'-',end5:4,end3:1,mismatches:1,thermo_viable:false}],binding_sites_truncated:50})`);
assert.ok(sites.includes('<details open>') && sites.includes('− ←'));
assert.ok(sites.includes('&lt;chr1>') && sites.includes('tools.no') && sites.includes('tools.omittedSites'));
const grouped = run(`renderCheck({primers:{F:'ACGA'},input_sequence_forms:{F:{input_5to3:'ACGA',reverse:'AGCA',complement_3to5:'TGCT',reverse_complement_5to3:'TCGT'}},results:[{oligos:{F:'TCGT'},reverse_complemented_inputs:['F'],products:[]},{oligos:{F:'ACGA'},reverse_complemented_inputs:[],products:[]}]})`);
assert.ok(grouped.indexOf('<h2>tools.literalGroup') < grouped.indexOf('<h2>tools.alternativeGroup'));
for (const sequence of ['ACGA','AGCA','TGCT','TCGT']) assert.ok(grouped.includes(sequence));
assert.ok(grouped.includes('tools.complement') && grouped.includes('tools.reverseComplement'));
const gated = run(`thermoCheckEvidence({thermo_status:'evaluated_defaults_gated',thermo_site_stats:{gated_per_primer:{F:2,R:1}}})`);
assert.ok(gated.includes('tools.thermoGated') && gated.includes('tools.thermoRejected: 3'));
const skipped = run(`thermoCheckEvidence({thermo_status:'skipped_no_associated_genome'})`);
assert.ok(skipped.includes('tools.thermoSkipped'));
const noSites = run(`thermoCheckEvidence({thermo_status:'evaluated_defaults_gated',thermo_evaluated:false})`);
assert.ok(noSites.includes('tools.thermoNoSites') && !noSites.includes('tools.thermoGated'));
assert.ok(run('I18N.ja["tools.oligo"]').startsWith('OFF'));
assert.ok(run('I18N.ja["tools.auto"]').startsWith('ON'));
assert.ok(run('I18N.ja["tools.orientationHint"]').includes('F/R・R/F・F/F・R/R'));
console.log('Standalone BLAST/Primer3, orientation hypotheses and FASTA evidence passed.');
