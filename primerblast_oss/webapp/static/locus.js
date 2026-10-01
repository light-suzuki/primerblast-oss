/* Coordinate-aware views. Coordinates come from the engine, never sequence search. */
Object.assign(I18N.ja, {
  'map.title': 'どこを増幅する？', 'map.template': '入力配列', 'map.target': '目的領域',
  'map.product': '増幅範囲', 'map.binding': '配列上の結合位置', 'map.reference': '参照配列（テンプレート方向）',
  'map.forward': 'Fプライマー →', 'map.reverse': '← Rプライマー',
  'map.coordinates': '表示座標は1始まり・両端を含みます。矢印は伸長方向です。',
  'map.noSequence': 'この保存結果には参照配列がありません。再計算すると結合位置を表示できます。',
  'map.offtarget': 'オフターゲットは何 bp？', 'map.gap': '標的産物とのサイズ差',
  'map.noOfftarget': '記録されたオフターゲット産物はありません。探索の完了状態も確認してください。',
  'map.engineered': '設計オリゴと参照配列の違いを色で表示します。dCAPSの意図的な変更を含みます。',
  'map.digest': '制限酵素でどこが切れる？', 'map.cut': '上鎖の切断位置（産物先頭からのbp）',
  'map.uncut': '切断なし',
  'map.changes': '変更塩基（テンプレート方向・参照 → 設計）',
  'map.gelOnly': '標的の切断パターンだけの模式図です。背景産物・実測泳動は含みません。特異性の未完了状態は変わりません。',
  'map.tiny': '破線は50 bp未満の断片です。検出や分離を保証しません。',
  'map.gel': 'ゲル電気泳動の模式図',
  'save.title': 'データの保存・再読込み', 'save.local': 'このブラウザーに保存',
  'save.export': 'プロジェクトJSONを出力', 'save.import': 'プロジェクトJSONを読込',
  'save.hint': '入力条件と結果を一緒に保存します。ブラウザー内の保存は直近5件。JSONは別の端末でも読み込めます。',
  'save.empty': '保存した結果を選択', 'save.delete': '選択した保存を削除',
  'save.success': '入力条件と結果をこのブラウザーに保存しました。',
  'save.failed': '保存・読込みができませんでした。容量やJSONの形式を確認してください。',
  'save.restored': '保存した結果を表示しています。計算は実行していません。',
  'save.needsResult': '計算が完了すると保存できます。', 'save.list': '保存したプロジェクト',
  'save.json': '保存用JSON', 'save.download': 'JSONファイルをダウンロード',
  'save.copyHint': 'ダウンロードできないブラウザーでは、下のJSONをコピーして .json ファイルとして保存できます。'
});
Object.assign(I18N.en, {
  'map.title': 'What will amplify?', 'map.template': 'Input sequence', 'map.target': 'Target region',
  'map.product': 'Amplicon', 'map.binding': 'Primer binding positions', 'map.reference': 'Reference sequence (template direction)',
  'map.forward': 'Forward primer →', 'map.reverse': '← Reverse primer',
  'map.coordinates': 'Displayed coordinates are 1-based, inclusive. Arrows indicate extension direction.',
  'map.noSequence': 'This saved result lacks reference sequence. Recalculate to show binding positions.',
  'map.offtarget': 'How long are the off-target products?', 'map.gap': 'Size difference from target',
  'map.noOfftarget': 'No off-target products were recorded. Check whether the search is complete.',
  'map.engineered': 'Differences between the designed oligo and reference are highlighted, including intentional dCAPS changes.',
  'map.digest': 'Where does the enzyme cut?', 'map.cut': 'Top-strand cut boundaries (bp from amplicon start)',
  'map.uncut': 'Uncut',
  'map.changes': 'Changed bases (template direction, reference → designed)',
  'map.gelOnly': 'Illustration of target digestion only. Background products and measured migration are excluded. Unresolved specificity remains unresolved.',
  'map.tiny': 'Dashed bands are below 50 bp. Detection and separation are not guaranteed.',
  'map.gel': 'Illustrative gel electrophoresis',
  'save.title': 'Save and restore your data', 'save.local': 'Save in this browser',
  'save.export': 'Export project JSON', 'save.import': 'Load project JSON',
  'save.hint': 'Keep inputs and results together. Browser storage retains five recent saves. JSON can be opened on another device.',
  'save.empty': 'Choose a saved result', 'save.delete': 'Delete selected save',
  'save.success': 'Inputs and results were saved in this browser.',
  'save.failed': 'Could not save or load. Check storage capacity and the JSON format.',
  'save.restored': 'Showing a saved result. No calculation was run.',
  'save.needsResult': 'Saving becomes available after calculation finishes.', 'save.list': 'Saved projects',
  'save.json': 'Project JSON', 'save.download': 'Download JSON file',
  'save.copyHint': 'If your browser cannot download, copy the JSON below into a .json file.'
});
function localSpan(span, context, genomic) {
  if (!span || span.length !== 2) return null;
  const sign = context.strand === '-' ? -1 : 1;
  const values = span.map(value => genomic ? (Number(value) - Number(context.anchor)) * sign : Number(value));
  return values.every(Number.isFinite) ? [Math.min(...values), Math.max(...values)] : null;
}
function displayCoordinate(index, context) {
  return Number(context.anchor || 1) + index * (context.strand === '-' ? -1 : 1);
}
function complement(sequence) {
  const codes = {A:'T',T:'A',C:'G',G:'C',R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',V:'B',D:'H',H:'D',N:'N'};
  return sequence.toUpperCase().split('').reverse().map(base => codes[base] || 'N').join('');
}
function bindingWindow(sequence, span, oligo, reverse, context) {
  const start = Math.max(0, span[0] - 12), end = Math.min(sequence.length - 1, span[1] + 12);
  const expected = reverse ? complement(oligo) : oligo.toUpperCase();
  const changes = [];
  let bases = '';
  for (let index = start; index <= end; index++) {
    const bound = index >= span[0] && index <= span[1];
    const changed = bound && expected[index - span[0]] && sequence[index].toUpperCase() !== expected[index - span[0]];
    if (changed) changes.push(`${displayCoordinate(index, context)}: ${sequence[index].toUpperCase()} → ${expected[index - span[0]]}`);
    bases += `<span class="${bound ? reverse ? 'base-r' : 'base-f' : ''}${changed ? ' base-changed' : ''}">${esc(sequence[index])}</span>`;
  }
  return `<div class="binding-window"><strong>${esc(t(reverse ? 'map.reverse' : 'map.forward'))}</strong><p>${esc(displayCoordinate(span[0], context))}–${esc(displayCoordinate(span[1], context))}</p><code>${esc(displayCoordinate(start, context))} ${bases} ${esc(displayCoordinate(end, context))}</code><p class="hint">5′-${esc(oligo)}-3′</p>${changes.length ? `<p>${esc(t('map.changes'))}: ${esc(changes.join(', '))}</p>` : ''}</div>`;
}
function locusView(context, forwardSpan, reverseSpan, pair, targetSpan) {
  const sequence = context.sequence || '';
  const length = sequence.length || Number(context.length);
  if (!length || !forwardSpan || !reverseSpan) return `<p class="hint">${esc(t('map.noSequence'))}</p>`;
  const all = [...forwardSpan, ...reverseSpan];
  const low = Math.min(...all), high = Math.max(...all);
  if (low < 0 || high >= length) return `<p class="hint">${esc(t('map.noSequence'))}</p>`;
  const x = index => 80 + 730 * index / Math.max(1, length - 1);
  const rect = (span, y, color) => `<rect x="${x(span[0])}" y="${y}" width="${Math.max(2, x(span[1]) - x(span[0]))}" height="13" rx="3" fill="${color}"/>`;
  const arrow = (span, y, color, reverse) => {
    const left = x(span[0]), right = x(span[1]);
    const tip = reverse ? left : right, tail = reverse ? right : left;
    return `<line x1="${tail}" y1="${y}" x2="${tip}" y2="${y}" stroke="${color}" stroke-width="7"/><path d="M${tip + (reverse ? 8 : -8)},${y - 8} L${tip},${y} L${tip + (reverse ? 8 : -8)},${y + 8}" fill="none" stroke="${color}" stroke-width="3"/>`;
  };
  let svg = `<svg class="locus-map" viewBox="0 0 900 190" role="img" aria-label="${esc(t('map.title'))}"><title>${esc(t('map.title'))}</title><line x1="80" y1="36" x2="810" y2="36" stroke="#bcc9c1" stroke-width="5"/><text x="80" y="20">${esc(displayCoordinate(0, context))}</text><text x="810" y="20" text-anchor="end">${esc(displayCoordinate(length - 1, context))}</text>`;
  if (targetSpan) svg += rect(targetSpan, 29, '#c4a65c');
  svg += rect([low, high], 82, '#1b6255') + `<text x="80" y="73">${esc(t('map.product'))} · ${esc(pair.product_size || high - low + 1)} bp</text>`;
  svg += arrow(forwardSpan, 130, '#1b6255', false) + arrow(reverseSpan, 160, '#315b91', true);
  svg += `<text x="${x(forwardSpan[0])}" y="115">F →</text><text x="${x(reverseSpan[1])}" y="183" text-anchor="end">← R</text></svg>`;
  const label = context.chrom ? `${context.chrom} · ${context.strand || '+'}` : t('map.template');
  let html = `<section class="locus-card"><h4>${esc(t('map.title'))} · ${esc(label)}</h4>${svg}<p class="hint">${esc(t('map.coordinates'))}</p><p>${esc(t('map.product'))}: ${esc(displayCoordinate(low, context))}–${esc(displayCoordinate(high, context))} · ${esc(pair.product_size || high - low + 1)} bp</p>`;
  if (targetSpan) html += `<p><span class="target-key"></span>${esc(t('map.target'))}: ${esc(displayCoordinate(targetSpan[0], context))}–${esc(displayCoordinate(targetSpan[1], context))}</p>`;
  if (typeof geneView === 'function') html += geneView(context, [low, high]);
  if (sequence) {
    const product = sequence.slice(low, high + 1);
    const name = `${context.chrom || 'template'}:${displayCoordinate(low, context)}-${displayCoordinate(high, context)}_reference_${context.strand || '+'}`;
    const fasta = `>${name}\n${product.match(/.{1,80}/g).join('\n')}\n`;
    html += dl('fasta', 'amplified_region.fa', fasta);
  }
  if (sequence && typeof sequenceLevelView === 'function') html += sequenceLevelView(context, forwardSpan, reverseSpan, pair);
  if (sequence) html += `<details class="adv"><summary>${esc(t('map.binding'))}</summary><p class="hint">${esc(t('map.reference'))} · ${esc(t('map.engineered'))}</p>${bindingWindow(sequence, forwardSpan, pair.forward, false, context)}${bindingWindow(sequence, reverseSpan, pair.reverse, true, context)}</details>`;
  return html + '</section>';
}
function productSizes(perDb, targetSize) {
  let rows = '';
  for (const db of perDb || []) for (const product of db.off_target || (db.products || []).filter(product => product.on_target === false)) {
    rows += `<tr><td>${esc(db.db && db.db.split('/').pop())}</td><td>${esc(product.subject)}:${esc(product.start)}–${esc(product.end)}</td><td><strong>${esc(product.size)} bp</strong></td><td>${Math.abs(Number(product.size) - Number(targetSize))} bp</td></tr>`;
  }
  return `<section class="offtarget-sizes"><h4>${esc(t('map.offtarget'))}</h4>${rows ? `<div class="tbl-wrap"><table><thead><tr><th>${esc(t('res.databases'))}</th><th>${esc(t('res.subject'))}</th><th>bp</th><th>${esc(t('map.gap'))}</th></tr></thead><tbody>${rows}</tbody></table></div>` : `<p class="hint">${esc(t('map.noOfftarget'))}</p>`}</section>`;
}
function digestView(result, size) {
  if (!result || !size) return '';
  let svg = `<svg class="locus-map" viewBox="0 0 900 125" role="img" aria-label="${esc(t('map.digest'))}"><title>${esc(t('map.digest'))}</title>`;
  let labels = '';
  for (const [index, allele] of ['a','b'].entries()) {
    const cuts = (result[`allele_${allele}_cuts`] || []).filter(cut => cut.complete && cut.top_cut > 0 && cut.top_cut < size);
    const y = 35 + index * 60;
    svg += `<text x="15" y="${y + 5}">${index ? 'BB' : 'AA'}</text><line x1="80" x2="810" y1="${y}" y2="${y}" stroke="#1b6255" stroke-width="5"/>`;
    for (const cut of cuts) { const x = 80 + 730 * cut.top_cut / size; svg += `<line x1="${x}" x2="${x}" y1="${y - 12}" y2="${y + 12}" stroke="#b04b16" stroke-width="3"/>`; }
    labels += `<p>${index ? 'BB' : 'AA'}: ${esc(cuts.length ? cuts.map(cut => cut.top_cut).join(', ') : t('map.uncut'))} · ${esc((result[`allele_${allele}_fragments`] || []).join(' / '))} bp</p>`;
  }
  return `<section class="locus-card"><h4>${esc(t('map.digest'))} · ${esc(result.enzyme)} (${esc(result.recognition)})</h4>${svg}</svg><p class="hint">${esc(t('map.cut'))}</p>${labels}</section>`;
}
function candidateGelView(analysis) {
  if (!analysis || !analysis.lanes) return '';
  const sizes = Object.values(analysis.lanes).flat().map(band => Number(band.size)).filter(size => size > 0);
  if (!sizes.length) return '';
  const low = Math.log10(Math.min(...sizes)), high = Math.log10(Math.max(...sizes));
  const y = size => 65 + 330 * (high - Math.log10(size)) / Math.max(.1, high - low);
  let svg = `<svg class="locus-map" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 900 440" role="img" aria-label="${esc(t('map.gel'))}"><title>${esc(t('map.gelOnly'))}</title><rect x="70" y="45" width="780" height="375" rx="6" fill="#18262b"/>`;
  for (const [index, name] of ['M','AA','AB','BB'].entries()) {
    const x = 135 + index * 190;
    svg += `<text x="${x + 40}" y="28" text-anchor="middle" fill="#172e25">${name}</text>`;
    for (const band of analysis.lanes[name] || []) {
      if (!(Number(band.size) > 0)) continue;
      const tiny = band.size < 50;
      svg += `<line x1="${x}" x2="${x + 80}" y1="${y(band.size)}" y2="${y(band.size)}" stroke="#e9f8ef" stroke-width="4"${tiny ? ' stroke-dasharray="8 5" opacity=".55"' : ''}/>`;
      if (name === 'M') svg += `<text x="${x + 85}" y="${y(band.size) + 5}" fill="#e9f8ef" style="fill:#e9f8ef;font-size:12px">${esc(band.size)}</text>`;
    }
  }
  svg += '</svg>';
  const labels = ['AA','AB','BB'].map(name => `<p>${name}: ${esc((analysis.lanes[name] || []).map(band => band.size).join(' / '))} bp</p>`).join('');
  return `<section class="locus-card"><h4>${esc(t('map.gel'))}</h4><p class="evidence-note">${esc(t('map.gelOnly'))}</p>${svg}${labels}<p class="hint">${esc(t('map.tiny'))} · ${esc(analysis.ladder || '')} · ${esc(analysis.gel_percent)}%</p>${dl('svg', 'candidate_digest_illustration.svg', svg)}</section>`;
}
const designWithoutMap = renderDesign;
renderDesign = function(data) {
  let html = designWithoutMap(data);
  for (const template of data.templates || []) for (const [index, pair] of (template.pairs || []).entries()) {
    const context = template.template || {sequence: template.template_sequence || '', length: template.template_len, anchor: 1, strand: '+'};
    const forward = [pair.left_start, pair.left_start + pair.forward.length - 1];
    const reverse = [pair.right_start - pair.reverse.length + 1, pair.right_start];
    html += `<details class="pair-map"${index === 0 ? ' open' : ''}><summary>${esc(template.template_id)} · #${index + 1} · ${esc(t('map.title'))}</summary>${locusView(context, forward, reverse, pair)}${productSizes(pair.specificity && pair.specificity.per_db, pair.product_size)}</details>`;
  }
  return html;
};
const sequenceWithoutMap = renderSequence;
renderSequence = function(data) {
  let html = sequenceWithoutMap(data);
  for (const plan of data.plans || [data]) for (const pair of plan.amplicons || []) {
    const context = plan.template || {length: Math.max(...(plan.region || [0,0])) + 1, anchor: 1};
    html += `<details class="pair-map" open><summary>${esc(plan.template_id)} · #${esc(pair.index)} · ${esc(t('map.title'))}</summary>${locusView(context, pair.forward_pos, pair.reverse_pos, pair, plan.region)}${productSizes(pair.specificity && pair.specificity.per_db, pair.product_size)}</details>`;
  }
  return html;
};
const assayWithoutMap = renderAssay;
renderAssay = function(data) {
  let html = assayWithoutMap(data);
  const context = data.template || {};
  const genomic = Boolean(context.chrom);
  const target = data.target && localSpan([data.target.start, data.target.end], context, genomic);
  for (const [index, pair] of (data.pairs || []).entries()) {
    const forward = localSpan(pair.forward_pos, context, genomic);
    const reverse = localSpan(pair.reverse_pos, context, genomic);
    html += `<details class="pair-map"${index === 0 ? ' open' : ''}><summary>#${index + 1} · ${esc(t('map.title'))}</summary>${locusView(context, forward, reverse, pair, target)}${productSizes(pair.per_db_products, pair.product_size)}</details>`;
    const derived = pair.caps && pair.caps.dcaps && pair.caps.dcaps.best;
    if (pair.caps) html += digestView(derived ? derived.digest : pair.caps.best_result, derived ? derived.product_size : pair.product_size);
    if ((!data.exports || !data.exports.svg) && pair.caps) {
      const digest = derived ? derived.digest : pair.caps.best_result;
      html += candidateGelView(digest && digest.gel_analysis);
    }
    if (typeof restrictionViews === 'function' && pair.caps) html += restrictionViews(pair, context);
    if (derived && derived.specificity) {
      const checked = derived.specificity;
      html += `<details class="pair-map" open><summary>#${index + 1} · dCAPS · ${esc(derived.enzyme)} (${esc(derived.recognition)})</summary>${locusView(context, localSpan(checked.forward_pos, context, genomic), localSpan(checked.reverse_pos, context, genomic), derived, target)}${productSizes(checked.per_db_products, derived.product_size)}</details>`;
    }
  }
  return html;
};

let projectParams = null;
let projectExportUrl = null;
function projectMessage(key) { const node = $('#project-message'); node.dataset.i18n = key; node.textContent = t(key); }
function captureProjectInput(mode, form) { projectParams = {mode, params: buildParams(mode, form)}; }
const runWithoutProject = runMode;
runMode = async function(mode, form) {
  await runWithoutProject(mode, form);
  syncProjectButtons();
};
function projectBundle() {
  if (!lastResult) throw new Error(t('save.needsResult'));
  return {schema: 'primerblast-oss-project', version: 1, saved_at: new Date().toISOString(), input: projectParams, result: lastResult};
}
function projectList() {
  try { const list = JSON.parse(localStorage.getItem('pb_projects') || '[]'); return Array.isArray(list) ? list : []; }
  catch (_) { return []; }
}
function syncProjectButtons() {
  ['#project-save','#project-export'].forEach(selector => { const button = $(selector); if (button) button.disabled = !lastResult || Boolean(studioJob); });
}
function refreshProjects() {
  const select = $('#project-list'); if (!select) return;
  select.setAttribute('aria-label', t('save.list'));
  select.innerHTML = `<option value="">${esc(t('save.empty'))}</option>`;
  projectList().forEach((project, index) => {
    const option = document.createElement('option'); option.value = index;
    option.textContent = `${project.saved_at} · ${project.result && project.result.mode || ''}`; select.append(option);
  });
}
function restoreProject(project) {
  if (studioJob) throw new Error(t('btn.running'));
  if (!project || project.schema !== 'primerblast-oss-project' || project.version !== 1 || !project.result || !['design','check','tile','sequence','assay','markers','makedb','blast','primer3'].includes(project.result.mode)) throw new Error(t('save.failed'));
  if (!project.input || project.input.mode !== project.result.mode || !project.input.params || typeof project.input.params !== 'object') throw new Error(t('save.failed'));
  const previous = lastResult;
  try { renderResult(project.result); }
  catch (error) { if (previous) renderResult(previous); throw error; }
  projectParams = project.input; lastResult = project.result; studioError = null;
  if (project.input && project.input.params) {
    const form = $(`[data-mode="${project.input.mode}"]`);
    if (form) {
      for (const [name, value] of Object.entries(project.input.params)) {
        const field = form.elements.namedItem(name); if (!field) continue;
        if (field.type === 'checkbox') field.checked = Boolean(value); else field.value = Array.isArray(value) ? value.join('\n') : value;
      }
      ['source','target_kind'].forEach(name => { const field = form.elements.namedItem(name); if (field) field.dispatchEvent(new Event('change')); });
      const params = project.input.params; selectedDbs = Array.isArray(params.db) ? params.db : [];
      customDbs = [...new Set([...customDbs, ...selectedDbs])]; $('#db-genomes').value = params.db_genomes || ''; loadDatabases();
    }
  }
  showTab(project.result.mode); renderResult(lastResult); syncProjectButtons();
  projectMessage('save.restored');
}
document.addEventListener('DOMContentLoaded', () => {
  const controls = document.createElement('details'); controls.className = 'project-store adv'; controls.open = true;
  controls.innerHTML = `<summary data-i18n="save.title"></summary><p class="hint" data-i18n="save.hint"></p><div class="project-actions"><button type="button" id="project-save" class="ghost" data-i18n="save.local"></button><button type="button" id="project-export" class="ghost" data-i18n="save.export"></button><label class="ghost"><span data-i18n="save.import"></span><input id="project-import" type="file" accept=".json"></label><select id="project-list" aria-label="Saved projects"></select><button type="button" class="ghost" id="project-delete" data-i18n="save.delete"></button></div><p id="project-message" role="status"></p>`;
  $('#results').insertAdjacentElement('beforebegin', controls);
  const output = document.createElement('details'); output.id = 'project-json-output'; output.hidden = true;
  output.innerHTML = `<summary data-i18n="save.json"></summary><p class="hint" data-i18n="save.copyHint"></p><a id="project-download" class="ghost" download="primerblast-project.json" data-i18n="save.download"></a><textarea id="project-json" readonly aria-label="Project JSON" rows="5"></textarea>`;
  controls.append(output);
  $('#project-save').addEventListener('click', () => {
    try { const list = projectList(); list.unshift(projectBundle()); localStorage.setItem('pb_projects', JSON.stringify(list.slice(0,5))); refreshProjects(); projectMessage('save.success'); }
    catch (_) { projectMessage('save.failed'); }
  });
  $('#project-export').addEventListener('click', () => {
    try {
    if (studioJob) throw new Error();
    const json = JSON.stringify(projectBundle(), null, 2);
    if (projectExportUrl) URL.revokeObjectURL(projectExportUrl);
    projectExportUrl = URL.createObjectURL(new Blob([json], {type:'application/json'}));
    $('#project-json').textContent = json; $('#project-json').value = json; output.hidden = false; output.open = true;
    const link = $('#project-download'); link.href = projectExportUrl; link.click();
    } catch (_) { projectMessage('save.failed'); }
  });
  $('#project-import').addEventListener('change', async event => {
    try { const file = event.target.files[0]; if (!file) return; if (file.size > 20 * 1024 * 1024) throw new Error(); restoreProject(JSON.parse(await file.text())); }
    catch (_) { projectMessage('save.failed'); }
    event.target.value = '';
  });
  $('#project-list').addEventListener('change', event => { if (event.target.value !== '') try { restoreProject(projectList()[Number(event.target.value)]); } catch (_) { projectMessage('save.failed'); } });
  $('#project-delete').addEventListener('click', () => { try { const index = $('#project-list').value; if (index === '') return; const list = projectList(); list.splice(Number(index),1); localStorage.setItem('pb_projects',JSON.stringify(list)); refreshProjects(); } catch (_) { projectMessage('save.failed'); } });
  const languageHook = onLangChange;
  onLangChange = function() { languageHook(); refreshProjects(); syncProjectButtons(); };
  applyI18n(); refreshProjects(); syncProjectButtons();
});
