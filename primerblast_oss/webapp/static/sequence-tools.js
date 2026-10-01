/* Standalone local tools and reference product exports. */
Object.assign(I18N.ja, {
  'tools.blast': '通常BLAST', 'tools.primer3': 'Primer3',
  'tools.blastHint': 'DNAまたは複数配列のFASTAを入力します。両鎖を検索し、相同性・座標・アラインメントを表示します。',
  'tools.primer3Hint': '配列だけでプライマーを設計します。検索DBは不要です。特異性は未確認なので、候補をPCR確認へ送って確認できます。',
  'tools.task': '検索方法', 'tools.limit': 'クエリごとの最大ヒット配列数',
  'tools.orientation': '入力配列の種類', 'tools.oligo': '発注済みオリゴ（5′→3′・両鎖を自動検索）',
  'tools.auto': 'ゲノムからコピーした配列（逆相補も自動確認・2本まで）',
  'tools.orientationHint': 'F/Rの入れ替えや＋/−鎖の指定は不要です。発注済み配列はそのまま入力。コピー配列では各入力と逆相補の組合せを別々に計算し、候補オリゴを表示します。単純な逆順や相補だけには変換しません。',
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
  'tools.blast': 'BLAST', 'tools.primer3': 'Primer3',
  'tools.blastHint': 'Enter DNA or multi-record FASTA. Search both strands and inspect identity, coordinates and alignments.',
  'tools.primer3Hint': 'Design primers from sequence alone, without a database. Specificity is not evaluated; send candidates to PCR check.',
  'tools.task': 'Search task', 'tools.limit': 'Maximum target sequences per query',
  'tools.orientation': 'Input sequence type', 'tools.oligo': 'Ordered oligos (5′→3′; both strands searched)',
  'tools.auto': 'Copied genomic sequence (also try reverse complements; up to two primers)',
  'tools.orientationHint': 'No F/R swap or +/− strand setting is needed. Paste ordered oligos unchanged. For copied sequences, each original/reverse-complement combination is calculated separately and candidate oligos are shown. No plain reverse or complement-only conversion is applied.',
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
renderCheck = function(data) {
  let html = dl('json', 'pcr_check_all.json', JSON.stringify(data, null, 2));
  for (const result of data.results || []) {
    html += checkWithoutSequence({primers: result.oligos || data.primers, results: [result]});
    html += `<p class="hint">${esc(t('studio.search'))}: ${esc(result.search_completeness || 'unknown')}</p>`;
    if (data.input_orientation === 'auto') html += `<p class="evidence-note">${esc(t('tools.hypothesis'))}: ${esc((result.reverse_complemented_inputs || []).join(', ') || t('tools.original'))}</p>`;
    if (result.fasta) html += dl('fasta', 'predicted_products.fa', result.fasta);
    for (const product of result.products || []) {
      html += `<details class="pair-card"><summary>${esc(product.subject)}:${product.start}–${product.end} · ${product.size} bp · ${esc(product.orientation)}</summary>`;
      html += primerRows({[product.fwd_primer]: (result.oligos || data.primers)[product.fwd_primer], [product.rev_primer]: (result.oligos || data.primers)[product.rev_primer]});
      const context = {chrom: product.subject, start: product.start, end: product.end, anchor: product.start, strand: '+', length: product.size, sequence: product.sequence || '', annotations: product.annotations};
      html += geneView(context, [0, product.size - 1]);
      if (product.sequence) html += `<h4>${esc(t('tools.reference'))}</h4><p class="hint">${esc(t('tools.referenceHint'))}</p>${dl('fasta', 'product.fa', product.fasta)}<textarea readonly rows="5" aria-label="FASTA">${esc(product.fasta)}</textarea>`;
      else html += `<p class="hint">${esc(t('tools.missing'))}</p><pre class="ascii">${esc(product.sequence_error || '')}</pre>`;
      html += '</details>';
    }
  }
  return html || checkWithoutSequence(data);
};

document.addEventListener('DOMContentLoaded', () => {
  const previousShowTab = showTab;
  showTab = name => { previousShowTab(name); if (name === 'primer3') $('.db-panel').hidden = true; };
  document.addEventListener('click', event => {
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
