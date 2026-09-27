const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const context = vm.createContext({
  I18N: {ja:{}, en:{}}, document:{addEventListener() {}},
  t: key => key, esc: value => String(value).replaceAll('<','&lt;')
});
vm.runInContext(fs.readFileSync('primerblast_oss/webapp/static/enzymes.js','utf8'), context);
const info = {pattern:{top:'5′ G|AATT C 3′',bottom:'3′ C TTAA|G 5′'},cuts:[[1,5]]};
const html = context.enzymePatternHTML(info);
assert.ok(html.includes('enz.over5'));
assert.ok(html.includes('4 nt'));
assert.ok(html.includes('3′ C TTAA|G 5′'));
assert.ok(context.enzymePatternHTML({...info,cuts:[[5,1]]}).includes('enz.over3'));
assert.ok(context.enzymePatternHTML({...info,cuts:[[3,3]]}).includes('enz.blunt'));
assert.ok(context.enzymePatternHTML({}).includes('enz.unknown'));
const names = context.enzymeRelationshipsHTML({same_cut_enzymes:['HpaII'],different_cut_enzymes:['<different>']});
assert.ok(names.includes('enz.same'));
assert.ok(names.includes('enz.different'));
assert.ok(names.includes('&lt;different>'));
assert.ok(!names.includes('<different>'));
console.log('Both-strand patterns, end types and enzyme relationship display passed.');
