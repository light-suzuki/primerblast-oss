/* Run with: node tests/web_locus.test.js (no npm dependencies). */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const context = vm.createContext({
  I18N: {ja: {}, en: {}}, document: {addEventListener() {}},
  renderDesign() { return ''; }, renderSequence() { return ''; },
  renderAssay() { return ''; }, runMode() {},
  t: key => key, esc: value => String(value ?? '').replaceAll('<', '&lt;'),
  dl() { return ''; },
});
vm.runInContext(fs.readFileSync('primerblast_oss/webapp/static/locus.js', 'utf8'), context);
vm.runInContext(fs.readFileSync('primerblast_oss/webapp/static/sequence-view.js', 'utf8'), context);
const run = expression => vm.runInContext(expression, context);
assert.equal(run('JSON.stringify(localSpan([180,199], {anchor:219,strand:"-"}, true))'), '[20,39]');
assert.equal(run('JSON.stringify(localSpan([100,119], {anchor:80,strand:"+"}, true))'), '[20,39]');
assert.equal(run('displayCoordinate(20, {anchor:219,strand:"-"})'), 199);
assert.equal(run('complement("AGCTRY")'), 'RYAGCT');
assert.equal(run('localSpan(["bad",119], {anchor:80},true)'), null);
const sizes = run('productSizes([{db:"ref",products:[{on_target:true,size:100},{on_target:false,subject:"chr2",start:20,end:141,size:122}]}],100)');
assert.ok(sizes.includes('122 bp'));
assert.ok(sizes.includes('22 bp'));
assert.ok(!sizes.includes('undefined'));
assert.ok(!sizes.includes('100 bp'));
assert.ok(run('productSizes([{db:"ref",off_target:[{subject:"chr3",start:1,end:300,size:300}]}],100)').includes('200 bp'));
const binding = run('bindingWindow("AAACCC", [0,2], "TTA", false, {anchor:1})');
assert.equal((binding.match(/base-changed/g) || []).length, 2);
const reverse = run('bindingWindow("AAACCC", [3,5], "GGG", true, {anchor:1})');
assert.equal((reverse.match(/base-changed/g) || []).length, 0);
const cuts = run('digestView({enzyme:"EcoRI",recognition:"GAATTC",allele_a_cuts:[{complete:true,top_cut:100}],allele_b_cuts:[],allele_a_fragments:[100,200],allele_b_fragments:[300]},300)');
assert.ok(cuts.includes('AA: 100'));
assert.ok(cuts.includes('BB: map.uncut'));
assert.ok(cuts.includes('100 / 200 bp'));
const gel = run('candidateGelView({lanes:{AA:[{size:214}],AB:[{size:214},{size:195},{size:19}],BB:[{size:195},{size:19}],M:[{size:500}]},ladder:"1kb",gel_percent:4})');
assert.ok(gel.includes('map.gelOnly'));
assert.ok(gel.includes('stroke-dasharray'));
assert.ok(gel.includes('BB: 195 / 19 bp'));
assert.equal(run('overlapBases([10,29],[20,39])'), 10);
const cleavage = run('cleavageWindow("GGGAATTCCC", {enzyme:"EcoRI",recognition:"GAATTC",site_pos:2,site_strand:"+",top_cut:3,bottom_cut:7,complete:true})');
assert.equal((cleavage.match(/cut-boundary/g)||[]).length, 2);
assert.ok(cleavage.includes('GAATTC (GAATTC, +)'));
assert.ok(cleavage.includes('↑ 3 / ↓ 7'));
const uncut = run('cleavageWindow("GGGACTTCCC", {enzyme:"EcoRI",recognition:"GAATTC",site_pos:2,site_strand:"+",top_cut:3,bottom_cut:7,complete:true},false)');
assert.ok(!uncut.includes('cut-boundary'));
assert.ok(uncut.includes('seq.compare'));
assert.ok(run('geneView({annotations:{status:"seqid_not_found"}},[0,10])').includes('seq.badSeqid'));
const page = run('sequencePage({context:{sequence:"A".repeat(801),anchor:1000,strand:"-"},forward:[10,29],reverse:[110,129]},600)');
assert.ok(page.includes('601–801 / 801 bp'));
assert.ok(page.includes('400'));
const gene = run('geneView({sequence:"A".repeat(140),anchor:219,strand:"-",annotations:{status:"loaded",source:"ref.gff3",genes:[{id:"g1",name:"G1",start:100,end:199,strand:"-",transcripts:[{id:"tx1",segments:[{type:"exon",start:100,end:120},{type:"exon",start:150,end:199}]}]}]}},[20,49])');
assert.ok(gene.includes('seq.exons: 1'));
assert.ok(gene.includes('seq.overlap: 30 bp'));
const engineeredPage = run('sequencePage({context:{sequence:"AACCGG",anchor:1,strand:"+"},forward:[0,1],reverse:[4,5],pair:{forward:"TA",reverse:"CC"}},0)');
assert.equal((engineeredPage.match(/base-changed/g)||[]).length, 1);
assert.ok(run('cleavageWindow("AAAAGGTCTCAAAAA", {enzyme:"BsaI",recognition:"GGTCTC",site_pos:4,site_strand:"+",top_cut:11,bottom_cut:15,complete:true})').includes('↑ 11 / ↓ 15'));
console.log('Binding orientation, coordinate conversion, off-target sizes and digest geometry passed.');
