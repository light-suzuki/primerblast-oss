/* Standalone local tools and reference product exports. */
Object.assign(I18N.ja, {
  'tools.showFasta':'FASTAを表示・コピー', 'tools.copyFasta':'FASTAをコピー',
  'tools.copied':'FASTAをコピーしました。', 'tools.copyFallback':'配列を選択しました。Ctrl+C / Command+Cでコピーしてください。',
  'tools.literalGroup': '入力配列そのものの検索結果',
  'tools.alternativeGroup': '入力を変更した配列の候補（別オリゴ）',
  'tools.forms': '入力配列と変換配列を見比べる',
  'tools.reverse': '逆順（文字の並びを反転）',
  'tools.complement': '相補配列（入力に対応する3′→5′）',
  'tools.reverseComplement': '逆相補配列（5′→3′）',
  'tools.formsHint': '入力したオリゴは自動で書き換えません。変換配列は照合用です。自動確認ONで検索する変更候補は逆相補です。逆順・相補配列は照合用に表示します。',
  'tools.reviewDirection': '配列の向きを確認：入力のままの2本の組合せでは産物が未検出ですが、逆相補へ変更した2本では候補があります。コピー元や発注配列と見比べてください。',
  'tools.blast': '通常BLAST', 'tools.primer3': 'Primer3',
  'tools.blastHint': 'DNAまたは複数配列のFASTAを入力します。両鎖を検索し、相同性・座標・アラインメントを表示します。',
  'tools.primer3Hint': '配列だけでプライマーを設計します。検索DBは不要です。特異性は未確認なので、候補をPCR確認へ送って確認できます。',
  'tools.task': '検索方法', 'tools.limit': 'クエリごとの最大ヒット配列数',
  'tools.orientation': '逆相補候補の自動確認', 'tools.oligo': 'OFF：入力配列のまま（両鎖を検索）',
  'tools.auto': 'ON：入力と逆相補の候補を比較（2本まで）',
  'tools.orientationHint': '配列は5′→3′で入力します。OFFは入力したオリゴだけを評価し、ONは逆相補に変更した別オリゴの候補も比較します。どちらも両鎖を検索し、F/R・R/F・F/F・R/Rの増幅候補と実際の伸長方向を表示します。入力配列は書き換えません。',
  'tools.genome': '増幅配列を取得するゲノムFASTA（任意・最初の検索DBと同じ参照）',
  'tools.gff3': '増幅領域の遺伝子を表示するGFF3（任意・最初の検索DBと同じ参照）',
  'tools.reference': '増幅領域の参照配列（＋鎖・1始まり両端を含む）',
  'tools.referenceHint': 'ゲノム由来の配列です。プライマーのミスマッチや5′テールを反映した実際のPCR産物配列ではありません。',
  'tools.missing': '配列を取得できません。対応するFASTAと索引（.fai）を指定してください。',
  'tools.hypothesis': '逆相補を使う候補入力', 'tools.original': '入力した配列のまま',
  'tools.check': 'このオリゴをPCR確認へ', 'tools.unscreened': '特異性は未確認',
  'tools.limited': 'ヒット数が上限に達しています。上限を増やして再検索できます。',
  'tools.hspHint': 'BLASTが返したアラインメントです。ヒットがない場合も相同性の不存在を証明するものではありません。',
  'dl.fasta': '参照FASTA', 'dl.tsv': 'BLAST TSV'
});
Object.assign(I18N.ja, {'f.num_threads': '検索に使うスレッド数', 'tools.annotationMissing': '遺伝子注釈を取得できませんでした。遺伝子の有無は未確認です。', 'tools.query': '1. 検索したい配列を入力', 'tools.searchDb': '2. 相同性を検索するゲノムを選択'});
Object.assign(I18N.en, {
  'tools.showFasta':'Show / copy FASTA', 'tools.copyFasta':'Copy FASTA',
  'tools.copied':'FASTA copied.', 'tools.copyFallback':'Sequence selected. Press Ctrl+C / Command+C to copy.',
  'tools.literalGroup': 'Search results for the supplied sequences',
  'tools.alternativeGroup': 'Candidates using changed sequences (different oligos)',
  'tools.forms': 'Compare input and transformed sequences',
  'tools.reverse': 'Reverse (character order reversed)',
  'tools.complement': 'Complement (3′→5′, aligned with the input)',
  'tools.reverseComplement': 'Reverse complement (5′→3′)',
  'tools.formsHint': 'Supplied oligos are never overwritten. These forms are for comparison. ON searches reverse-complement alternatives. Reverse and complement-only forms are displayed for comparison.',
  'tools.reviewDirection': 'Review sequence direction: no two-primer product was observed as supplied, but reverse-complement changes yield candidates. Compare the source or ordered oligos.',
  'tools.blast': 'BLAST', 'tools.primer3': 'Primer3',
  'tools.blastHint': 'Enter DNA or multi-record FASTA. Search both strands and inspect identity, coordinates and alignments.',
  'tools.primer3Hint': 'Design primers from sequence alone, without a database. Specificity is not evaluated; send candidates to PCR check.',
  'tools.task': 'Search task', 'tools.limit': 'Maximum target sequences per query',
  'tools.orientation': 'Automatic reverse-complement alternatives', 'tools.oligo': 'OFF: supplied oligos only (both strands searched)',
  'tools.auto': 'ON: compare original and reverse-complement alternatives (up to two primers)',
  'tools.orientationHint': 'Enter sequences 5′→3′. OFF evaluates the supplied oligos; ON also compares separately changed reverse-complement oligos. Both modes search both strands and all F/R, R/F, F/F and R/R combinations, with actual extension directions. Supplied sequences are never overwritten.',
  'tools.genome': 'Genome FASTA for product extraction (optional; matches the first search DB)',
  'tools.gff3': 'GFF3 for overlapping genes (optional; matches the first search DB)',
  'tools.reference': 'Reference product sequence (+ strand; 1-based inclusive)',
  'tools.referenceHint': 'Genomic reference bases, without oligo mismatches or 5′ tails incorporated into the actual PCR product.',
  'tools.missing': 'Sequence unavailable. Supply the matching FASTA and its .fai index.',
  'tools.hypothesis': 'Inputs reverse-complemented for this candidate', 'tools.original': 'As supplied',
  'tools.check': 'Send these oligos to PCR check', 'tools.unscreened': 'Specificity not evaluated',
  'tools.limited': 'Target limit reached. Increase the limit to search again.',
  'tools.hspHint': 'Alignments reported by BLAST. No hits does not establish the absence of homology.',
  'dl.fasta': 'Reference FASTA', 'dl.tsv': 'BLAST TSV'
});
Object.assign(I18N.en, {'f.num_threads': 'Search threads', 'tools.annotationMissing': 'Annotation could not be loaded. Gene presence is unresolved.', 'tools.query': '1. Enter sequences to search', 'tools.searchDb': '2. Select genomes for homology search'});

Object.assign(I18N.ja, {
  'tools.inputVerdict': '入力配列のままでの増幅予測',
  'tools.literalProduct': '入力のままで産物を予測',
  'tools.changedProduct': '逆相補に変更した別オリゴの候補',
  'tools.noLiteral': '指定した検索条件では産物が見つかりませんでした。',
  'tools.unknownLiteral': '産物は未検出ですが、検索が不完全なため増幅可否は未確定です。',
  'tools.pairProducts': '異なる2本での産物', 'tools.selfProducts': '同じ1本による産物',
  'tools.alternatives': '逆相補へ変更した候補産物',
  'tools.errorUnknown': '入力ミスか、ゲノム表示からコピーした配列かは検索だけでは断定できません。入力と発注済みの5′→3′配列を照合してください。これは計算上の予測です。',
  'tools.swapped': 'F/Rラベルが逆（左がR、右がF）。ラベルの入れ替え自体は配列変更を意味しません。',
  'tools.self': '同じプライマーが両側に結合する産物です。F/Rペアの成功とは別です。',
  'tools.direction': '参照上の向き・伸長方向：＋ → 右、− ← 左。両端から内向きに伸びる組合せです。',
  'tools.inputSequence': '入力（5′→3′）', 'tools.usedSequence': 'この候補のオリゴ（5′→3′）',
  'tools.noOriginalLocus': 'この領域の産物は入力配列のままでは未検出です（指定した検索条件内）。',
  'tools.unknownOriginalLocus': 'この領域は入力のままでは未検出ですが、元の検索が不完全なため未確定です。',
  'tools.originalLocus': 'この領域は入力配列のままでも産物を検出しています。',
  'tools.sites': '各プライマーの一致候補・方向',
  'tools.siteHint': '一致候補だけでは増幅を意味しません。向き・距離・ミスマッチ・熱力学条件も判定します。座標は1始まり。＋/−は参照へのアラインメント方向で、結合する鋳型鎖は反対側です。',
  'tools.omittedSites': '省略した一致候補数', 'tools.mismatch': 'ミスマッチ数',
  'tools.thermoStatus': '熱力学評価',
  'tools.thermoSkipped': '対応するゲノムFASTAがないため未評価',
  'tools.thermoUnavailable': '計算ライブラリが利用できないため未評価',
  'tools.thermoDisabled': '評価を無効にしています',
  'tools.thermoNoSites': '一致候補がないため評価する部位がありません',
  'tools.thermoGated': '評価済み。条件に適合しない部位は増幅予測から除外',
  'tools.thermoAnnotation': '評価済み。部位の除外には使っていません',
  'tools.thermoPartial': '一部の部位を評価できませんでした。増幅可否は未確定です',
  'tools.thermoFailed': '部位を評価できませんでした。増幅可否は未確定です',
  'tools.thermoRejected': '熱力学条件で除外した部位数',
  'tools.thermo': '熱力学的な伸長可否', 'tools.yes': '適合', 'tools.no': '不適合', 'tools.unknown': '未評価'
});
Object.assign(I18N.en, {
  'tools.inputVerdict': 'Amplification prediction for the supplied sequences',
  'tools.literalProduct': 'Product predicted as supplied',
  'tools.changedProduct': 'Alternative oligo with reverse-complement changes',
  'tools.noLiteral': 'No products found under the selected search conditions.',
  'tools.unknownLiteral': 'No products observed; incomplete search leaves amplification unresolved.',
  'tools.pairProducts': 'Products using two different primers', 'tools.selfProducts': 'Products using one primer on both sides',
  'tools.alternatives': 'Alternative products after reverse-complement changes',
  'tools.errorUnknown': 'Search alone cannot establish a typo or a copied genomic sequence. Compare the input with the ordered 5′→3′ oligos. These are computational predictions.',
  'tools.swapped': 'F/R labels are swapped (R on the left, F on the right). Swapping labels does not itself change an oligo sequence.',
  'tools.self': 'The same primer binds both ends. This is separate from a successful F/R pair.',
  'tools.direction': 'Reference alignment / extension: + → right, − ← left. Both ends extend inward.',
  'tools.inputSequence': 'Input (5′→3′)', 'tools.usedSequence': 'Oligo for this candidate (5′→3′)',
  'tools.noOriginalLocus': 'No product at this locus was observed with the supplied oligos under the selected conditions.',
  'tools.unknownOriginalLocus': 'No product at this locus was observed as supplied; incomplete original search leaves this unresolved.',
  'tools.originalLocus': 'A product at this locus was also observed with the supplied oligos.',
  'tools.sites': 'Primer alignment candidates and directions',
  'tools.siteHint': 'An alignment alone does not imply amplification. Direction, distance, mismatches and thermodynamic conditions also matter. Coordinates are 1-based. +/− denotes reference alignment direction; the physical template strand is opposite.',
  'tools.omittedSites': 'Omitted alignment candidates', 'tools.mismatch': 'Mismatches',
  'tools.thermoStatus': 'Thermodynamic evaluation',
  'tools.thermoSkipped': 'Not evaluated: no associated genome FASTA',
  'tools.thermoUnavailable': 'Not evaluated: calculation library unavailable',
  'tools.thermoDisabled': 'Evaluation disabled',
  'tools.thermoNoSites': 'No candidate sites to evaluate',
  'tools.thermoGated': 'Evaluated; nonviable sites excluded from product predictions',
  'tools.thermoAnnotation': 'Evaluated for annotation only; sites are not excluded',
  'tools.thermoPartial': 'Some sites could not be evaluated; amplification remains unresolved',
  'tools.thermoFailed': 'Sites could not be evaluated; amplification remains unresolved',
  'tools.thermoRejected': 'Sites excluded by thermodynamic conditions',
  'tools.thermo': 'Thermodynamic viability', 'tools.yes': 'Viable', 'tools.no': 'Not viable', 'tools.unknown': 'Not evaluated'
});

function inputAssessments(data) {
  let html = '';
  for (const assessment of data.input_assessments || []) {
    html += `<section class="run-summary"><h3>${esc(t('tools.inputVerdict'))} · ${esc(assessment.db)}</h3>`;
    if (assessment.as_supplied_products) html += `<p class="evidence-note">${esc(t('tools.literalProduct'))}: ${assessment.as_supplied_products}</p>`;
    else html += `<p class="evidence-note">${esc(t(assessment.original_search_complete ? 'tools.noLiteral' : 'tools.unknownLiteral'))}</p>`;
    html += `<p>${esc(t('tools.pairProducts'))}: ${assessment.distinct_primer_products} · ${esc(t('tools.selfProducts'))}: ${assessment.same_primer_products} · ${esc(t('tools.alternatives'))}: ${assessment.reverse_complement_candidates}</p></section>`;
    if (!assessment.distinct_primer_products && assessment.distinct_primer_reverse_complement_candidates) html += `<p class="evidence-note">${esc(t('tools.reviewDirection'))}</p>`;
  }
  return html + `<p class="hint">${esc(t('tools.errorUnknown'))}</p>`;
}

function inputSequenceForms(data) {
  if (!Object.keys(data.input_sequence_forms || {}).length) return '';
  let html = `<details open><summary>${esc(t('tools.forms'))}</summary><p class="hint">${esc(t('tools.formsHint'))}</p>`;
  for (const [name, forms] of Object.entries(data.input_sequence_forms || {})) {
    html += `<h4>${esc(name)}</h4><div class="table-wrap"><table><tbody>`;
    for (const [key, label] of [['input_5to3', 'tools.inputSequence'], ['reverse', 'tools.reverse'], ['complement_3to5', 'tools.complement'], ['reverse_complement_5to3', 'tools.reverseComplement']]) html += `<tr><th>${esc(t(label))}</th><td class="seqmono">${esc(forms[key])}</td></tr>`;
    html += '</tbody></table></div>';
  }
  return html + '</details>';
}

function thermoCheckEvidence(result) {
  if (!result.thermo_status) return '';
  let label = {
    skipped_no_associated_genome: 'tools.thermoSkipped', unavailable: 'tools.thermoUnavailable',
    disabled: 'tools.thermoDisabled', evaluated_defaults_gated: 'tools.thermoGated',
    evaluated_gated: 'tools.thermoGated', evaluated_defaults_annotation_only: 'tools.thermoAnnotation',
    evaluated_annotation_only: 'tools.thermoAnnotation', partial_unresolved_sites: 'tools.thermoPartial',
    failed_no_resolvable_sites: 'tools.thermoFailed'
  }[result.thermo_status];
  if (result.thermo_evaluated === false && result.thermo_status.startsWith('evaluated')) label = 'tools.thermoNoSites';
  const rejected = Object.values((result.thermo_site_stats || {}).gated_per_primer || {}).reduce((a, b) => a + b, 0);
  return `<p class="hint">${esc(t('tools.thermoStatus'))}: ${esc(label ? t(label) : result.thermo_status)}${rejected ? ` · ${esc(t('tools.thermoRejected'))}: ${rejected}` : ''}</p>`;
}

function bindingSiteEvidence(result) {
  if (!result.binding_site_counts || !Object.keys(result.binding_site_counts).length) return '';
  let html = `<details${(result.products || []).length ? '' : ' open'}><summary>${esc(t('tools.sites'))}</summary><p class="hint">${esc(t('tools.siteHint'))}</p>`;
  for (const [name, counts] of Object.entries(result.binding_site_counts)) html += `<p>${esc(name)}: ＋ → ${counts['+']} · − ← ${counts['-']}</p>`;
  html += `<div class="table-wrap"><table><thead><tr><th>Primer</th><th>Subject</th><th>5′ → 3′</th><th>＋/−</th><th>${esc(t('tools.mismatch'))}</th><th>${esc(t('tools.thermo'))}</th></tr></thead><tbody>`;
  for (const site of result.binding_sites || []) html += `<tr><td>${esc(site.primer)}</td><td>${esc(site.subject)}</td><td>${site.end5} → ${site.end3}</td><td>${site.reference_strand === '+' ? '＋ →' : '− ←'}</td><td>${site.mismatches}</td><td>${esc(t(site.thermo_viable === true ? 'tools.yes' : site.thermo_viable === false ? 'tools.no' : 'tools.unknown'))}</td></tr>`;
  html += '</tbody></table></div>';
  if (result.binding_sites_truncated) html += `<p>${esc(t('tools.omittedSites'))}: ${result.binding_sites_truncated}</p>`;
  return html + '</details>';
}

function productInputEvidence(product) {
  const evidence = product.input_evidence;
  if (!evidence) return '';
  let html = `<p class="evidence-note">${esc(t(evidence.sequence_status === 'as_supplied' ? 'tools.literalProduct' : 'tools.changedProduct'))}</p>`;
  if (evidence.changed_primers.length) html += `<p>${esc(t(evidence.as_supplied_locus_observed ? 'tools.originalLocus' : evidence.original_search_complete ? 'tools.noOriginalLocus' : 'tools.unknownOriginalLocus'))}</p>`;
  if (evidence.label_orientation === 'swapped') html += `<p>${esc(t('tools.swapped'))}</p>`;
  if (evidence.label_orientation === 'same_primer') html += `<p>${esc(t('tools.self'))}</p>`;
  html += `<p class="hint">${esc(t('tools.direction'))}</p><pre class="ascii">${esc(product.fwd_primer)} ＋ →     ${product.start} … ${product.end}     ← − ${esc(product.rev_primer)}</pre>`;
  for (const binding of product.primer_bindings || []) {
    html += `<p><strong>${esc(binding.primer)} ${binding.reference_strand === '+' ? '＋ →' : '− ←'}</strong> · 5′ ${binding.end5} / 3′ ${binding.end3 == null ? '?' : binding.end3} · ${esc(t('tools.mismatch'))}: ${binding.mismatches == null ? '?' : binding.mismatches}</p>`;
    html += `<p>${esc(t('tools.inputSequence'))}: <span class="seqmono">${esc(binding.input_sequence)}</span></p>`;
    if (binding.oligo !== binding.input_sequence) html += `<p class="evidence-note">${esc(t('tools.usedSequence'))}: <span class="seqmono">${esc(binding.oligo)}</span></p>`;
  }
  return html;
}

function renderBlast(data) {
  let html = `<h3>${esc(t('tools.blast'))} · ${esc(data.task)}</h3><p class="hint">${esc(t('tools.hspHint'))}</p>` + dl('json', 'blast.json', JSON.stringify(data, null, 2));
  const names = Object.fromEntries((data.queries || []).map(q => [q.id, q.name]));
  for (const result of data.results || []) {
    html += `<section class="run-summary"><h3>${esc(result.db)}</h3>${dl('tsv', 'blast.tsv', result.tsv)}`;
    if ((result.at_target_limit || []).length) html += `<p class="evidence-note">${esc(t('tools.limited'))} ${esc(result.at_target_limit.map(id => names[id] || id).join(', '))}</p>`;
    if (!result.hits.length) html += `<p>${esc(t('res.none'))}</p>`;
    for (const hit of result.hits) html += `<details class="pair-card"><summary>${esc(names[hit.qseqid] || hit.qseqid)} → ${esc(hit.sseqid)}:${hit.sstart}–${hit.send} (${esc(hit.strand)}) · ${hit.pident}% · ${hit.length} bp · E=${hit.evalue}</summary><p>Query: ${hit.qstart}–${hit.qend} · bit score ${hit.bitscore}</p><pre class="ascii">Query   ${esc(hit.qseq)}\nSubject ${esc(hit.sseq)}</pre></details>`;
    html += '</section>';
  }
  return html;
}

let primer3Candidates = [];
function renderPrimer3(data) {
  primer3Candidates = [];
  let html = `<h3>Primer3</h3><p class="evidence-note">${esc(t('tools.unscreened'))}</p>` + dl('json', 'primer3.json', JSON.stringify(data, null, 2));
  for (const template of data.templates || []) {
    html += `<section class="run-summary"><h3>${esc(template.template_id)}</h3>`;
    for (const pair of template.pairs || []) {
      const index = primer3Candidates.push(pair) - 1;
      html += `<details class="pair-card" open><summary>#${pair.index + 1} · ${pair.product_size} bp · Tm ${pair.tm_f.toFixed(1)} / ${pair.tm_r.toFixed(1)} °C</summary>${primerRows({F: pair.forward, R: pair.reverse})}<p>GC ${pair.gc_f.toFixed(1)} / ${pair.gc_r.toFixed(1)}% · penalty ${pair.penalty.toFixed(2)}</p><button class="ghost" type="button" data-primer3-check="${index}">${esc(t('tools.check'))}</button>${dl('fasta', 'primer3_product.fa', pair.fasta)}${locusView({sequence: template.template_sequence, anchor: 1, strand: '+'}, [pair.left_start, pair.left_start + pair.left_len - 1], [pair.right_start - pair.right_len + 1, pair.right_start], pair)}</details>`;
    }
    if (!template.pairs.length) html += `<p>${esc(t('res.none'))}</p>`;
    html += `<details><summary>${esc(t('res.explain'))}</summary><pre class="ascii">${esc(template.primer3_explain)}</pre></details></section>`;
  }
  return html;
}

const checkWithoutSequence = renderCheck;
function fastaCopyPanel(fasta) {
  return `<details class="fasta-export"><summary>${esc(t('tools.showFasta'))}</summary><div class="fasta-copy"><button class="ghost" type="button" data-copy-fasta>${esc(t('tools.copyFasta'))}</button><p role="status"></p><textarea readonly rows="5" aria-label="FASTA">${esc(fasta)}</textarea></div></details>`;
}
Object.assign(I18N.ja, {
  'tools.length':'予測断片長（bp）', 'tools.export':'参照FASTAを保存',
  'tools.noFragment':'予測産物が0件のため、断片長と産物FASTAはありません。一致候補だけから産物配列を作ることはできません。',
  'tools.extractFailed':'この産物の参照配列は未取得です。詳細に取得エラーを表示しています。'
});
Object.assign(I18N.en, {
  'tools.length':'Predicted fragment length (bp)', 'tools.export':'Save reference FASTA',
  'tools.noFragment':'No predicted products: there is no product length or FASTA to export. Binding candidates alone do not define a product sequence.',
  'tools.extractFailed':'Reference sequence is unavailable for this product. Extraction errors are shown in its details.'
});
renderCheck = function(data) {
  let html = (data.partial ? `<p class="evidence-note">${esc(t('run.partial'))} ${data.completed_units}/${data.total_units}</p>` : inputAssessments(data)) + inputSequenceForms(data) + dl('json', 'pcr_check_all.json', JSON.stringify(data, null, 2));
  const literal = (data.results || []).filter(r => !(r.reverse_complemented_inputs || []).length);
  const alternatives = (data.results || []).filter(r => (r.reverse_complemented_inputs || []).length);
  for (const [label, results] of [['tools.literalGroup', literal], ['tools.alternativeGroup', alternatives]]) {
    if (!results.length) continue;
    html += `<h2>${esc(t(label))}</h2>`;
    for (const result of results) {
      html += checkWithoutSequence({primers: result.oligos || data.primers, results: [result]});
      if (!(result.products || []).length) html += `<p class="evidence-note">${esc(t('tools.noFragment'))}</p>`;
      html += `<p class="hint">${esc(t('studio.search'))}: ${esc(result.search_completeness || 'unknown')}</p>`;
      html += thermoCheckEvidence(result);
      if ((result.reverse_complemented_inputs || []).length) html += `<p class="evidence-note">${esc(t('tools.hypothesis'))}: ${esc(result.reverse_complemented_inputs.join(', '))}</p>`;
      html += bindingSiteEvidence(result);
      if (result.fasta) html += dl('fasta', 'predicted_products.fa', result.fasta) + fastaCopyPanel(result.fasta);
      for (const product of result.products || []) {
        const status = product.input_evidence && t(product.input_evidence.sequence_status === 'as_supplied' ? 'tools.literalProduct' : 'tools.changedProduct');
        html += `<details class="pair-card"><summary>${esc(product.subject)}:${product.start}–${product.end} · ${product.size} bp · ${esc(product.fwd_primer)} ＋ → / ${esc(product.rev_primer)} − ←${status ? ' · ' + esc(status) : ''}</summary>`;
        html += productInputEvidence(product);
        html += primerRows({[product.fwd_primer]: (result.oligos || data.primers)[product.fwd_primer], [product.rev_primer]: (result.oligos || data.primers)[product.rev_primer]});
        const context = {chrom: product.subject, start: product.start, end: product.end, anchor: product.start, strand: '+', length: product.size, sequence: product.sequence || '', annotations: product.annotations};
        html += geneView(context, [0, product.size - 1]);
        if (product.sequence) html += `<h4>${esc(t('tools.reference'))}</h4><p class="hint">${esc(t('tools.referenceHint'))}</p>${dl('fasta', 'product.fa', product.fasta)}${fastaCopyPanel(product.fasta)}`;
        else html += `<p class="hint">${esc(t('tools.missing'))}</p><pre class="ascii">${esc(product.sequence_error || '')}</pre>`;
        html += '</details>';
      }
    }
  }
  return html || checkWithoutSequence(data);
};

document.addEventListener('DOMContentLoaded', () => {
  const previousShowTab = showTab;
  showTab = name => { previousShowTab(name); if (name === 'primer3') $('.db-panel').hidden = true; };
  document.addEventListener('click', event => {
    const copyButton = event.target.closest('[data-copy-fasta]');
    if (copyButton) {
      const panel = copyButton.closest('.fasta-copy');
      const textarea = panel.querySelector('textarea');
      const message = panel.querySelector('[role="status"]');
      Promise.resolve().then(() => navigator.clipboard.writeText(textarea.value))
        .then(() => { message.textContent = t('tools.copied'); })
        .catch(() => { textarea.focus(); textarea.select(); message.textContent = t('tools.copyFallback'); });
      return;
    }
    const button = event.target.closest('[data-primer3-check]');
    if (!button) return;
    const pair = primer3Candidates[Number(button.dataset.primer3Check)];
    if (!pair) return;
    const form = $('[data-mode="check"]');
    form.elements.forward.value = pair.forward;
    form.elements.reverse.value = pair.reverse;
    form.elements.primers.value = '';
    form.elements.input_orientation.value = 'as_supplied';
    showTab('check');
  });
  applyI18n();
});
