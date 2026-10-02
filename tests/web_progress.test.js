/* Progress updates preserve opened partial results and label measured stages. */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const nodes = Object.fromEntries(['#results','#job-progress','#job-partial'].map(key => [key,{innerHTML:'',insertAdjacentHTML(where,html){this.innerHTML += html;}}]));
const context = vm.createContext({
  I18N:{ja:{},en:{}}, document:{addEventListener(){}},
  renderResult(){},renderSequence(){},renderAssay(){},onLangChange(){},runMode(){},
  $:key=>nodes[key], $$:()=>[], _dlStore:[],
  t:key=>key, esc:value=>String(value ?? '').replaceAll('<','&lt;').replaceAll('>','&gt;'),
  dbNames:dbs=>dbs.map(value=>value.split('/').pop()).join(', '),
  renderCheck:data=>'partial:'+data.completed_units, wireDownloads(){},
  Date:{now:()=>10000}
});
vm.runInContext(fs.readFileSync('primerblast_oss/webapp/static/studio.js','utf8'),context);
const run=code=>vm.runInContext(code,context);
const timing=run(`runTiming({total_seconds:12.34,stage_seconds:{blast:2.5,realign:9.84}})`);
assert.ok(timing.includes('12.3') && timing.includes('run.blast') && timing.includes('9.8'));
run(`studioJob={start:0}; studioObserve({progress:{stage:'realign',database:'x/<ref>',elapsed_seconds:10,stage_elapsed_seconds:3,completed:249,total:1000,completed_units:1,total_units:4,hypothesis:2,hypotheses:4,primer:'F'},partial_revision:1,partial_result:{completed_units:1}})`);
assert.ok(nodes['#job-progress'].innerHTML.includes('249/1000'));
assert.ok(nodes['#job-progress'].innerHTML.includes('1/4'));
assert.ok(nodes['#job-progress'].innerHTML.includes('&lt;ref&gt;'));
assert.equal(nodes['#job-partial'].innerHTML,'partial:1');
nodes['#job-partial'].innerHTML='opened details retained';
run('studioRunning()');
assert.equal(nodes['#job-partial'].innerHTML,'opened details retained');
run(`studioObserve({progress:{stage:'thermo'},partial_revision:1,partial_result:{completed_units:1}})`);
assert.equal(nodes['#job-partial'].innerHTML,'opened details retained');
run(`studioObserve({progress:{stage:'export'},partial_revision:2,partial_result:{completed_units:2}})`);
assert.equal(nodes['#job-partial'].innerHTML,'partial:2');
run(`studioShowError({message:'second DB failed',partial:{completed_units:1}})`);
assert.ok(nodes['#results'].innerHTML.includes('studio.failed'));
assert.ok(nodes['#results'].innerHTML.includes('partial:1'));
console.log('Stage timing, partial revisions and open-detail preservation passed.');
