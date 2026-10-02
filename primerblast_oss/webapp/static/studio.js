/* Offline, bilingual interaction layer. Engine evidence remains unchanged. */
const studioWords = {
  ja: {
    'tab.sequence': '配列を読み切る', 'studio.skip': '作業画面へ移動',
    'studio.step1': '目的を選ぶ', 'studio.step2': '配列と検索先を指定', 'studio.step3': '候補と根拠を確認',
    'studio.fastaMap': '検索DBに対応するゲノムFASTA（任意）',
    'studio.fastaMapHint': '1行に DBパス=FASTAパス を指定します。パスはこのサーバー側のものです。対応FASTAがあればプライマー全長を再確認できます。',
    'studio.sequenceWhen': '長い領域を重複するPCR断片で読む',
    'studio.sequenceIntro': '目的領域を重複する増幅断片に分け、未カバー領域と発注用配列を確認します。PCRによる領域カバーは、Sangerの実測読取り長を保証しません。',
    'studio.strand': '対象の鎖方向', 'studio.ampliconSize': '増幅断片のサイズ範囲（bp）',
    'studio.m13': '発注用配列にM13テールを付ける', 'studio.pending': '判定保留',
    'studio.pendingHint': '探索が未完了、または全長の確認が不足しています。特異的とはまだ判断できません。',
    'studio.file': 'FASTA / テキストを読み込む', 'studio.fileHint': '配列を入力欄に読み込みます。ゲノムDBのアップロードではありません（最大2 MB）。',
    'studio.pathHint': 'ファイルパスはサーバーから読める場所を指定してください。WSLでは /mnt/c/... 形式です。',
    'studio.invalid': '入力を確認してください', 'studio.required': '必須項目が空です。',
    'studio.dna': 'DNA配列（ACGTとIUPAC文字）またはFASTAを入力してください。',
    'studio.interval': 'chr1:100-1000 の形式で、開始≦終了の正の座標を指定してください。',
    'studio.range': '小さい値-大きい値 の形式で正の範囲を指定してください。',
    'studio.number': '数値の範囲を確認してください。', 'studio.failed': '計算を完了できませんでした',
    'studio.retry': '入力、サーバー側のファイルパス、BLAST / Primer3の準備状況を確認してください。',
    'studio.details': '技術的な詳細', 'studio.elapsed': '経過', 'studio.seconds': '秒',
    'studio.computing': '計算中です。検索DBの大きさによって時間が変わります。',
    'studio.evidence': '計算上の候補です。実際のPCR・消化・配列読取りによる検証は別途必要です。',
    'studio.coverage': 'PCR領域カバー', 'studio.gaps': '未カバー区間', 'studio.complete': '全領域をカバー',
    'studio.amplicons': '増幅断片', 'studio.order': '発注用配列（5′→3′）',
    'studio.region': 'カバーする座標（0始まり・両端を含む）', 'studio.overlap': '前の断片との重複',
    'studio.marker': 'マーカー解析の根拠', 'studio.noEvidence': 'この候補には該当する解析結果がありません。',
    'studio.localFile': 'ファイルを読み込めませんでした。テキスト形式で2 MB以下のファイルを選んでください。'
  },
  en: {
    'tab.sequence': 'Sequence a region', 'studio.skip': 'Skip to workspace',
    'studio.step1': 'Choose your goal', 'studio.step2': 'Add sequence and databases', 'studio.step3': 'Review candidates and evidence',
    'studio.fastaMap': 'Genome FASTA for each search database (optional)',
    'studio.fastaMapHint': 'Enter one DB path=FASTA path per line. Paths belong to the server. A matching FASTA enables full-length primer verification.',
    'studio.sequenceWhen': 'Read a long region using overlapping PCR fragments',
    'studio.sequenceIntro': 'Split the target into overlapping amplicons, review coverage gaps, and export ordering sequences. PCR coverage does not guarantee measured Sanger read coverage.',
    'studio.strand': 'Target strand', 'studio.ampliconSize': 'Amplicon size range (bp)',
    'studio.m13': 'Add M13 tails to ordering sequences', 'studio.pending': 'Unresolved',
    'studio.pendingHint': 'The search is incomplete or full-length verification is missing. Specificity cannot yet be established.',
    'studio.file': 'Load FASTA / text', 'studio.fileHint': 'Loads sequence into the text box, rather than uploading a genome database (maximum 2 MB).',
    'studio.pathHint': 'Use paths accessible to the server. Under WSL, use /mnt/c/... paths.',
    'studio.invalid': 'Check your input', 'studio.required': 'This field is required.',
    'studio.dna': 'Enter DNA (ACGT and IUPAC characters) or FASTA.',
    'studio.interval': 'Use chr1:100-1000 with positive coordinates and start ≤ end.',
    'studio.range': 'Enter a positive range as minimum-maximum.',
    'studio.number': 'Check the numeric range.', 'studio.failed': 'The calculation could not finish',
    'studio.retry': 'Check inputs, server file paths, and BLAST / Primer3 availability.',
    'studio.details': 'Technical details', 'studio.elapsed': 'Elapsed', 'studio.seconds': 's',
    'studio.computing': 'Calculating. Runtime depends on database size.',
    'studio.evidence': 'These are computational candidates. PCR, digestion, and sequencing still require experimental validation.',
    'studio.coverage': 'PCR region coverage', 'studio.gaps': 'Uncovered intervals', 'studio.complete': 'Full region coverage',
    'studio.amplicons': 'Amplicons', 'studio.order': 'Ordering sequences (5′→3′)',
    'studio.region': 'Covered coordinates (0-based, inclusive)', 'studio.overlap': 'Overlap with previous fragment',
    'studio.marker': 'Marker analysis evidence', 'studio.noEvidence': 'No corresponding analysis is available for this candidate.',
    'studio.localFile': 'Could not read this file. Choose a text file no larger than 2 MB.'
  }
};
Object.keys(studioWords).forEach(lang => Object.assign(I18N[lang], studioWords[lang]));
Object.assign(I18N.ja, {'dl.svg': '仮想ゲル SVG', 'studio.bands': '予測バンド（bp）', 'studio.oligos': 'マーカー用プライマー（5′→3′）', 'studio.scorable': '計算上は判別可能', 'studio.ambiguous': '判別困難', 'studio.search': '検索・判定の状態'});
Object.assign(I18N.en, {'dl.svg': 'Virtual gel SVG', 'studio.bands': 'Predicted bands (bp)', 'studio.oligos': 'Marker primers (5′→3′)', 'studio.scorable': 'Computationally distinguishable', 'studio.ambiguous': 'Ambiguous', 'studio.search': 'Search and assessment status'});
Object.assign(I18N.ja, {'studio.screened': '探索完了（詳細で各プライマー組を確認）', 'studio.gel': 'ゲル・マーカー解析設定', 'studio.ladder': 'DNAラダー', 'studio.auto': '自動', 'studio.gelPercent': 'アガロース濃度（%）', 'studio.screenLimit': 'AS-PCR探索に使う親ペア数', 'studio.customBands': 'カスタムラダー（bp、カンマ区切り）'});
Object.assign(I18N.en, {'studio.screened': 'Search complete (review each primer set in details)', 'studio.gel': 'Gel and marker analysis settings', 'studio.ladder': 'DNA ladder', 'studio.auto': 'Automatic', 'studio.gelPercent': 'Agarose concentration (%)', 'studio.screenLimit': 'Parent pairs to screen for AS-PCR', 'studio.customBands': 'Custom ladder (bp, comma separated)'});
Object.assign(I18N.ja, {'studio.ampliconHint': '1つの範囲を指定します。例：500-800', 'studio.custom': 'カスタム'});
Object.assign(I18N.en, {'studio.ampliconHint': 'Enter one range, e.g. 500-800', 'studio.custom': 'Custom'});
Object.assign(I18N.ja, {'studio.overlapShort': '指定した重複長を下回ります。隣接する断片を確認してください。', 'studio.offline': 'サーバーに接続できません', 'studio.requestedOverlap': '指定した重複長'});
Object.assign(I18N.en, {'studio.overlapShort': 'Below the requested overlap. Review the neighboring amplicons.', 'studio.offline': 'Server unavailable', 'studio.requestedOverlap': 'Requested overlap'});
Object.assign(I18N.ja, {'studio.indexMissing': 'ゲノムFASTAの索引（.fai）が見つかりません。指定パスと索引の作成を確認してください。'});
Object.assign(I18N.en, {'studio.indexMissing': 'The genome FASTA index (.fai) is missing. Check the path and create the index.'});
let studioJob = null;
let studioError = null;
Object.assign(I18N.ja, {
  'run.prepare':'入力・参照の準備', 'run.blast':'BLAST検索', 'run.realign':'一致候補の参照配列での再確認',
  'run.thermo':'熱力学評価', 'run.enumerate':'産物の向き・距離の判定', 'run.export':'参照FASTA・遺伝子注釈の取得',
  'run.partial':'途中結果：評価が完了した検索単位のみ表示しています。残りは計算中です。',
  'run.units':'完了した検索単位', 'run.hypothesis':'入力候補', 'run.count':'この工程の処理件数',
  'run.total':'計算時間', 'run.breakdown':'工程別時間', 'run.phaseElapsed':'この工程の経過'
});
Object.assign(I18N.en, {
  'run.prepare':'Preparing input and references', 'run.blast':'Searching with BLAST', 'run.realign':'Rechecking candidates against reference sequence',
  'run.thermo':'Evaluating thermodynamic conditions', 'run.enumerate':'Checking product directions and distances', 'run.export':'Extracting lengths, FASTA and annotations',
  'run.partial':'Partial results: only completed search units are shown. Remaining units are still running.',
  'run.units':'Completed search units', 'run.hypothesis':'Input hypothesis', 'run.count':'Processed in this stage',
  'run.total':'Calculation time', 'run.breakdown':'Time by stage', 'run.phaseElapsed':'Elapsed in this stage'
});

function runTiming(timing) {
  if (!timing) return '';
  const seconds = value => Number(value || 0).toFixed(1);
  return `<section class="run-summary"><strong>${esc(t('run.total'))}: ${seconds(timing.total_seconds)} ${esc(t('studio.seconds'))}</strong><details><summary>${esc(t('run.breakdown'))}</summary>${Object.entries(timing.stage_seconds || {}).map(([stage, value]) => `<p>${esc(t('run.' + stage))}: ${seconds(value)} ${esc(t('studio.seconds'))}</p>`).join('')}</details></section>`;
}
const resultWithoutTiming = renderResult;
renderResult = function(data) {
  resultWithoutTiming(data);
  $('#results').insertAdjacentHTML('afterbegin', runTiming(data.timing));
};
function searchUnresolved(spec) {
  return !spec || typeof spec.specific_all_db !== 'boolean' || spec.search_complete_all_db === false ||
    (spec.search_completeness != null && spec.search_completeness !== 'complete');
}
function evidenceNote() { return `<p class="evidence-note">${esc(t('studio.evidence'))}</p>`; }
function renderSequence(data) {
  const plans = data.plans || (data.coverage ? [data] : []);
  let html = evidenceNote() + `<div class="downloads">${dl('json', 'sequence.json', JSON.stringify(data, null, 2))}${data.exports && data.exports.order ? dl('order', 'sequence_order.tsv', data.exports.order) : ''}</div>`;
  for (const plan of plans) {
    const c = plan.coverage || {};
    const pct = Number(c.coverage_fraction || 0) * 100;
    html += `<section class="run-summary"><h3>${esc(plan.template_id)}</h3><div class="metrics"><strong>${esc(t('studio.coverage'))}: ${pct.toFixed(1)}%</strong><span>${esc(t('studio.amplicons'))}: ${esc(plan.n_amplicons)}</span><span>${esc(t('studio.requestedOverlap'))}: ${esc(plan.requested_overlap)} bp</span></div><p>${esc(c.full_coverage ? t('studio.complete') : t('studio.gaps') + ': ' + JSON.stringify(c.gaps || []))}</p><div class="coverage-track"><span style="width:${Math.min(100, Math.max(0, pct))}%"></span></div>`;
    for (const a of plan.amplicons || []) {
      html += `<details class="pair-card" open><summary>#${esc(a.index)} · ${esc(a.product_size)} bp · ${esc(t('studio.overlap'))}: ${esc(a.overlap_to_prev)} bp</summary><p>${esc(t('studio.region'))}: ${esc((a.covers || []).join('–'))}</p><p>${esc(t('studio.order'))}</p><div class="primer-list"><code>F ${esc(a.order_forward || a.forward)}</code><code>R ${esc(a.order_reverse || a.reverse)}</code></div><p class="hint">${esc(t('studio.details'))}: ${esc(a.specificity && a.specificity.specificity_status || 'unknown')}</p></details>`;
      if (a.index > 1 && a.overlap_to_prev < plan.requested_overlap) html += `<p class="evidence-note">#${esc(a.index)} · ${esc(t('studio.overlapShort'))}</p>`;
      if (a.genomic) html += `<p class="hint">${esc(a.genomic.chrom)}:${esc(a.genomic.start)}–${esc(a.genomic.end)} (${esc(a.genomic.strand)})</p>`;
    }
    html += '</section>';
  }
  return html;
}
const originalAssayRenderer = renderAssay;
renderAssay = function(data) {
  let html = evidenceNote() + originalAssayRenderer(data);
  if (data.exports && data.exports.svg) html += `<div class="gel-preview"><img alt="${esc(t('dl.svg'))}" src="data:image/svg+xml;charset=utf-8,${encodeURIComponent(data.exports.svg)}">${dl('svg', 'virtual_gel.svg', data.exports.svg)}</div>`;
  for (const [index, pair] of (data.pairs || []).entries()) {
    const analyses = { CAPS: pair.caps, ASPCR: pair.aspcr };
    for (const [label, analysis] of Object.entries(analyses)) {
      if (!analysis) continue;
      html += `<section class="run-summary"><h3>#${index + 1} · ${esc(t('studio.marker'))} · ${label}</h3>`;
      if (label === 'CAPS') {
        const verdict = analysis.marker_verdict || '';
        const pending = !verdict || verdict.startsWith('indeterminate');
        const scorable = analysis.gel_scorable_all_db === true && !pending;
        html += `<p><span class="badge ${pending ? 'v-pending' : scorable ? 'v-ok' : 'v-bad'}">${esc(t(pending ? 'studio.pending' : scorable ? 'studio.scorable' : 'studio.ambiguous'))}</span> ${esc(analysis.best_marker_type || 'CAPS / dCAPS')} · ${esc(analysis.best_enzyme || '—')}</p>`;
        if (pending) html += `<p>${esc(t('studio.pendingHint'))}</p>`;
        html += `<p>${esc(t('studio.bands'))}</p><div class="marker-plot">`;
        for (const [name, fragments] of [['AA', analysis.allele_ref_fragments], ['BB', analysis.allele_alt_fragments]]) html += `<div class="marker-lane"><strong>${name}</strong><p>${esc((fragments || []).join(' / ') || '—')}</p></div>`;
        const lanes = analysis.gel_analysis && analysis.gel_analysis.lanes;
        if (lanes && lanes.AB) html += `<div class="marker-lane"><strong>AB</strong><p>${esc(lanes.AB.map(b => b.size).join(' / '))}</p></div>`;
        html += '</div>';
        const dcaps = analysis.dcaps && analysis.dcaps.best;
        if (dcaps) html += primerRows({F: dcaps.forward, R: dcaps.reverse});
      } else {
        const tetra = analysis.best_tetra;
        if (tetra) {
          html += `<h4>tetra-ARMS · ${esc(t('studio.oligos'))}</h4>` + primerRows({outer_F: tetra.outer_forward, outer_R: tetra.outer_reverse, ref_inner: tetra.ref_inner && tetra.ref_inner.primer, alt_inner: tetra.alt_inner && tetra.alt_inner.primer});
        }
        for (const allele of ['ref', 'alt']) {
          const candidate = analysis['best_classical_' + allele];
          if (candidate) html += `<h4>AS-PCR · ${allele} · ${esc(candidate.product_size)} bp</h4>` + primerRows({AS: candidate.primer, common: candidate.common_primer});
        }
        // A constructed oligo is not evidence of a completed specificity screen.
        const screen = analysis.specificity_screen || {};
        html += `<p class="evidence-note">${esc(t(screen.search_complete_all_sets === true ? 'studio.screened' : 'studio.pendingHint'))}</p>`;
      }
      html += `<details class="adv"><summary>${esc(t('studio.details'))}</summary><pre class="ascii">${esc(JSON.stringify(analysis, null, 2))}</pre></details></section>`;
    }
  }
  return html;
};
function primerRows(primers) {
  return `<div class="primer-list">${Object.entries(primers).filter(([, seq]) => seq).map(([name, seq]) => `<div class="primer-row"><strong>${esc(name)}</strong><code>${esc(seq)}</code></div>`).join('')}</div>`;
}
function studioRunning() {
  if (!$('#job-progress')) $('#results').innerHTML = '<div id="job-progress"></div><div id="job-partial"></div>';
  const p = studioJob.progress || {};
  const age = studioJob.observedAt ? (Date.now() - studioJob.observedAt) / 1000 : 0;
  const elapsed = p.elapsed_seconds == null ? (Date.now() - studioJob.start) / 1000 : p.elapsed_seconds + age;
  let html = `<div class="run-summary" role="status"><span class="spinner"></span><strong>${esc(t('btn.running'))}</strong><p>${esc(p.stage ? t('run.' + p.stage) : t('studio.computing'))}</p><p>${esc(t('studio.elapsed'))}: ${Math.floor(elapsed)} ${esc(t('studio.seconds'))}</p>`;
  if (p.database) html += `<p>${esc(dbNames([p.database]))} · ${esc(t('run.hypothesis'))}: ${p.hypothesis || 1}/${p.hypotheses || 1}${p.primer ? ' · ' + esc(p.primer) : ''}</p>`;
  if (p.total_units != null) html += `<p>${esc(t('run.units'))}: ${p.completed_units}/${p.total_units}</p>`;
  if (p.total != null) html += `<p>${esc(t('run.count'))}: ${p.completed || 0}/${p.total}</p>`;
  if (p.stage_elapsed_seconds != null) html += `<p>${esc(t('run.phaseElapsed'))}: ${Math.floor(p.stage_elapsed_seconds + age)} ${esc(t('studio.seconds'))}</p>`;
  $('#job-progress').innerHTML = html + '</div>';
}
function studioObserve(job) {
  studioJob.progress = job.progress;
  studioJob.observedAt = Date.now();
  studioRunning();
  if (job.partial_result && studioJob.partialRevision !== job.partial_revision) {
    studioJob.partialRevision = job.partial_revision;
    studioJob.partial = job.partial_result;
    _dlStore = [];
    $('#job-partial').innerHTML = renderCheck(job.partial_result);
    wireDownloads();
  }
}
function studioShowError(error) {
  const message = String(error.message || error.error || error);
  $('#results').innerHTML = `<div class="err"><strong>${esc(t('studio.failed'))}</strong><p>${esc(t(message.includes('FASTA index not found') ? 'studio.indexMissing' : 'studio.retry'))}</p><details><summary>${esc(t('studio.details'))}</summary><pre>${esc(message)}\n${esc(error.trace || '')}</pre></details></div>`;
  if (error.partial) {
    _dlStore = [];
    $('#results').insertAdjacentHTML('beforeend', renderCheck(error.partial));
    wireDownloads();
  }
}
function validateStudio(form, mode) {
  $$('.field-error', form).forEach(el => el.remove());
  $$('[aria-invalid]', form).forEach(el => { el.removeAttribute('aria-invalid'); el.removeAttribute('aria-describedby'); });
  const p = buildParams(mode, form);
  let first = null;
  function invalid(name, key) {
    const el = form.elements.namedItem(name);
    if (!el || el.getAttribute('aria-invalid') === 'true') return;
    el.setAttribute('aria-invalid', 'true');
    const note = document.createElement('span');
    note.className = 'field-error'; note.dataset.i18n = key; note.textContent = t(key);
    const id = `error-${mode}-${name}`; note.id = id; el.setAttribute('aria-describedby', id);
    el.insertAdjacentElement('afterend', note); first = first || el;
  }
  for (const [name, value] of Object.entries(p)) {
    const el = form.elements.namedItem(name);
    if (!el) continue;
    if (el.type === 'number' && (!Number.isFinite(Number(value)) || !el.checkValidity())) invalid(name, 'studio.number');
  }
  const required = mode === 'makedb' ? ['infile'] : [];
  if (mode === 'check' && !p.forward && !p.reverse && !(p.primers || []).length) invalid('forward', 'studio.required');
  if (mode === 'blast' || mode === 'primer3') required.push('template');
  if (mode === 'blast' && Number(p.evalue) <= 0) invalid('evalue', 'studio.number');
  if (mode === 'design' || mode === 'sequence') {
    const source = p.source || p.src || p.design_src || 'sequence';
    required.push(source === 'gene' ? 'gene' : source === 'interval' ? 'interval' : 'template');
    if (source === 'gene' || source === 'interval') required.push('genome');
    if (source === 'gene') required.push('gff3');
  }
  if (mode === 'assay') { required.push('genome', p.target_kind || 'gene'); if (p.target_kind === 'gene') required.push('gff3'); }
  if (mode === 'markers') required.push('genome', 'interval');
  if (mode === 'tile') required.push('template');
  for (const name of required) if (!p[name] || (Array.isArray(p[name]) && !p[name].length)) invalid(name, 'studio.required');
  if (p.template) {
    const dna = p.template.split(/\r?\n/).filter(line => !line.trim().startsWith('>')).join('').replace(/\s/g, '');
    if (!dna || /[^ACGTRYSWKMBDHVN]/i.test(dna)) invalid('template', 'studio.dna');
  }
  for (const name of ['forward', 'reverse']) if (p[name] && /[^ACGTRYSWKMBDHVN]/i.test(p[name].replace(/\s/g, ''))) invalid(name, 'studio.dna');
  if (p.interval) {
    const match = /^.+:(\d+)-(\d+)$/.exec(p.interval.trim());
    if (!match || +match[1] < 1 || +match[2] < +match[1]) invalid('interval', 'studio.interval');
  }
  for (const name of ['product_size', 'amplicon_size']) if (p[name]) {
    const ranges = name === 'product_size' ? p[name].split(',') : [p[name]];
    for (const range of ranges) {
      const match = /^(\d+)-(\d+)$/.exec(range.trim());
      if (!match || +match[1] < 1 || +match[2] < +match[1]) invalid(name, 'studio.range');
    }
  }
  if (first) { const details = first.closest('details'); if (details) details.open = true; first.focus(); }
  return !first;
}
runMode = async function(mode, form) {
  if (studioJob) return;
  if (!validateStudio(form, mode)) return;
  if (!['makedb', 'primer3'].includes(mode) && !selectedDbs.length) {
    $('#results').innerHTML = `<div class="err">${esc(t('db.pick'))}</div>`; $('#db-custom-input').focus(); return;
  }
  studioError = null; lastResult = null;
  if (typeof captureProjectInput === 'function') captureProjectInput(mode, form);
  studioJob = {start: Date.now(), mode};
  if (typeof syncProjectButtons === 'function') syncProjectButtons();
  $$('.run').forEach(btn => { btn.disabled = true; });
  $('#results').setAttribute('aria-busy', 'true'); studioRunning();
  const timer = setInterval(studioRunning, 1000);
  try {
    const response = await fetch(`/api/run/${mode}`, {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({params: buildParams(mode, form)})});
    const submission = await response.json();
    if (!response.ok || submission.error || !submission.job_id) throw new Error(submission.error || `HTTP ${response.status}`);
    let job;
    for (;;) {
      await new Promise(resolve => setTimeout(resolve, 1200));
      const poll = await fetch(`/api/job/${encodeURIComponent(submission.job_id)}`);
      if (!poll.ok) throw new Error(`HTTP ${poll.status}`);
      job = await poll.json();
      if (job.status === 'done' || job.status === 'error') break;
      studioObserve(job);
    }
    if (job.status === 'error') throw Object.assign(new Error(job.error), {trace: job.trace, timing: job.timing, partial: job.partial_result});
    lastResult = {...job.result, timing: job.timing};
  } catch (error) { error.partial = error.partial || studioJob.partial; studioError = error; }
  finally {
    clearInterval(timer); studioJob = null;
    $$('.run').forEach(btn => { btn.disabled = false; });
    $('#results').setAttribute('aria-busy', 'false');
    if (studioError) { studioShowError(studioError); $('#results').insertAdjacentHTML('afterbegin', runTiming(studioError.timing)); } else renderResult(lastResult);
  }
};
const originalLanguageHook = onLangChange;
onLangChange = function() {
  originalLanguageHook();
  document.documentElement.lang = LANG;
  document.title = LANG === 'ja' ? 'PrimerBLAST OSS — プライマー設計' : 'PrimerBLAST OSS — Primer design';
  if (studioJob) {
    studioRunning();
    if (studioJob.partial) {
      _dlStore = [];
      $('#job-partial').innerHTML = renderCheck(studioJob.partial);
      wireDownloads();
    }
  } else if (studioError) studioShowError(studioError);
};
document.addEventListener('DOMContentLoaded', () => {
  const tabs = $$('.tab');
  if (tabs.length) tabs[0].parentElement.setAttribute('role', 'tablist');
  tabs.forEach((tab, index) => {
    const mode = tab.dataset.tab;
    tab.id = `tab-${mode}`; tab.setAttribute('role', 'tab'); tab.setAttribute('aria-controls', `panel-${mode}`);
    const panel = $(`[data-panel="${mode}"]`);
    if (panel) { panel.id = `panel-${mode}`; panel.setAttribute('role', 'tabpanel'); panel.setAttribute('aria-labelledby', tab.id); }
    tab.addEventListener('keydown', event => {
      const visible = tabs.filter(el => getComputedStyle(el).display !== 'none');
      const current = visible.indexOf(tab);
      let next;
      if (event.key === 'ArrowRight') next = visible[(current + 1) % visible.length];
      if (event.key === 'ArrowLeft') next = visible[(current + visible.length - 1) % visible.length];
      if (event.key === 'Home') next = visible[0];
      if (event.key === 'End') next = visible[visible.length - 1];
      if (next) { event.preventDefault(); showTab(next.dataset.tab); next.focus(); }
    });
  });
  const source = $('#sequence-src');
  const sync = () => $$('[data-seq-src]').forEach(el => { el.style.display = el.dataset.seqSrc === source.value || (el.dataset.seqSrc === 'genomic' && source.value !== 'sequence') ? '' : 'none'; });
  source.addEventListener('change', sync); source.form.addEventListener('reset', () => setTimeout(sync, 0)); sync();
  $$('.panel[data-mode]').forEach(form => { form.noValidate = true; form.addEventListener('reset', () => {
    $$('.field-error', form).forEach(el => el.remove());
    $$('[aria-invalid]', form).forEach(el => { el.removeAttribute('aria-invalid'); el.removeAttribute('aria-describedby'); });
  }); });
  const assay = $('[data-mode="assay"]');
  const settings = document.createElement('details'); settings.className = 'adv';
  settings.innerHTML = `<summary data-i18n="studio.gel"></summary><fieldset><label><span data-i18n="studio.ladder"></span><select name="gel_ladder"><option value="auto" data-i18n="studio.auto"></option><option value="100bp">100 bp</option><option value="1kb">1 kb</option><option value="custom" data-i18n="studio.custom"></option></select></label><label><span data-i18n="studio.gelPercent"></span><select name="gel_percent"><option value="auto" data-i18n="studio.auto"></option><option>1</option><option>1.5</option><option>2</option><option>2.5</option><option>3</option><option>4</option></select></label><label><span data-i18n="studio.customBands"></span><input name="ladder_bands" placeholder="100,200,300,500,1000"></label><label><span data-i18n="studio.screenLimit"></span><input name="aspcr_pairs_to_screen" type="number" min="0" max="10" value="2"></label></fieldset>`;
  assay.querySelector('.actions').insertAdjacentElement('beforebegin', settings);
  $$('textarea[name="template"]').forEach(area => {
    const box = document.createElement('div'); box.className = 'file-tools';
    const group = area.closest('[data-src],[data-seq-src]');
    if (group && group.dataset.src) box.dataset.src = group.dataset.src;
    if (group && group.dataset.seqSrc) box.dataset.seqSrc = group.dataset.seqSrc;
    const label = document.createElement('label'); const caption = document.createElement('span'); caption.dataset.i18n = 'studio.file';
    const input = document.createElement('input'); input.type = 'file'; input.accept = '.fa,.fasta,.fna,.txt';
    label.append(caption, input); box.append(label);
    const hint = document.createElement('p'); hint.className = 'hint'; hint.dataset.i18n = 'studio.fileHint'; box.append(hint);
    area.closest('label').insertAdjacentElement('afterend', box);
    input.addEventListener('change', async () => {
      try { const file = input.files[0]; if (!file) return; if (file.size > 2 * 1024 * 1024) throw new Error(); area.value = await file.text(); }
      catch (_) { studioError = new Error(t('studio.localFile')); studioShowError(studioError); }
      input.value = '';
    });
  });
  $$('input[name="genome"],input[name="gff3"],input[name="fasta"]').forEach(input => {
    const hint = document.createElement('em'); hint.className = 'hint'; hint.dataset.i18n = 'studio.pathHint'; input.insertAdjacentElement('afterend', hint);
  });
  EXAMPLES.sequence = form => { setVal(form, 'template', DEMO_SEQ); setVal(form, 'template_id', 'example'); source.value = 'sequence'; sync(); };
  sync();
  applyI18n();
});
