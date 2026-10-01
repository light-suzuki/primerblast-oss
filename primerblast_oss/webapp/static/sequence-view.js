/* Sequence, annotation and cleavage evidence, using the engine's coordinates. */
Object.assign(I18N.ja, {
  'seq.title':'配列で増幅領域を確認', 'seq.start':'表示開始（配列内・1始まり）',
  'seq.show':'表示', 'seq.prev':'前へ', 'seq.next':'次へ',
  'seq.legend':'緑＋上線：F、青＋上線：R、灰色：増幅領域、橙＋太い下線：設計オリゴとの相違。細い下線は注釈上のエクソンです。',
  'seq.genes':'遺伝子のどこを増幅する？', 'seq.annotation':'表示用GFF3注釈（任意）',
  'seq.noAnnotation':'遺伝子注釈は未読込みです。ゲノム座標とGFF3を指定すると遺伝子・エクソン・CDSを表示できます。',
  'seq.noGenes':'読込み済みの注釈では、この表示領域に重なる遺伝子はありません。',
  'seq.badSeqid':'GFF3に一致する染色体名がありません。遺伝子がないという判定ではありません。',
  'seq.exons':'重なるエクソン', 'seq.overlap':'増幅領域との重なり',
  'seq.cutTitle':'認識配列と上下の鎖の切断位置', 'seq.enzymes':'切れる制限酵素（配列と切断部位）',
  'seq.cutLegend':'塗りつぶしは認識配列の位置（切れない側は比較位置）、│は塩基間の切断境界です。上下の鎖を同じ方向に揃えて表示します。座標はPCR産物内の1始まりです。',
  'seq.compare':'切断なし・比較用の同じ位置',
  'seq.noCut':'このアレルには記録された完全な切断部位がありません。',
});
Object.assign(I18N.en, {
  'seq.title':'Inspect the amplified sequence', 'seq.start':'Display start (1-based input index)',
  'seq.show':'Show', 'seq.prev':'Previous', 'seq.next':'Next',
  'seq.legend':'Green + top border: F; blue + top border: R; gray: amplicon; orange + thick underline: differences from the oligo. Thin underlines mark annotated exons.',
  'seq.genes':'Where in the gene does amplification occur?', 'seq.annotation':'GFF3 annotation for display (optional)',
  'seq.noAnnotation':'Gene annotations are not loaded. Supply genomic coordinates and GFF3 to show genes, exons and CDS.',
  'seq.noGenes':'The loaded annotation has no genes overlapping this displayed interval.',
  'seq.badSeqid':'No matching chromosome name in GFF3. This does not mean that genes are absent.',
  'seq.exons':'Overlapping exons', 'seq.overlap':'Overlap with amplicon',
  'seq.cutTitle':'Recognition sequence and cuts on both strands', 'seq.enzymes':'Enzymes that cut (sequence and cleavage sites)',
  'seq.cutLegend':'Shading marks recognition-site positions (comparison positions in the uncut allele); │ marks a boundary between bases. Both strands are aligned in the same direction. Coordinates are 1-based within the PCR product.',
  'seq.compare':'Uncut; matching position for comparison',
  'seq.noCut':'No complete cleavage site was recorded for this allele.',
});
let sequenceViews = [];
if (typeof renderResult === 'function') {
  const renderWithSequence = renderResult;
  renderResult = function(data) { sequenceViews = []; renderWithSequence(data); };
}
function overlapBases(a, b) { return Math.max(0, Math.min(a[1], b[1]) - Math.max(a[0], b[0]) + 1); }
function geneView(context, product) {
  const annotation = context.annotations;
  if (!annotation || annotation.status === 'not_provided' || annotation.status === 'no_genomic_coordinates') return `<p class="hint">${esc(t('seq.noAnnotation'))}</p>`;
  if (annotation.status === 'seqid_not_found') return `<p class="evidence-note">${esc(t('seq.badSeqid'))}</p>`;
  if (annotation.status === 'unavailable') return `<p class="evidence-note">${esc(t('tools.annotationMissing'))}</p>`;
  if (!(annotation.genes || []).length) return `<p class="hint">${esc(t('seq.noGenes'))}</p>`;
  const length = context.sequence.length || context.length || 0;
  let html = `<h4>${esc(t('seq.genes'))}</h4>`;
  for (const gene of annotation.genes) {
    const span = localSpan([gene.start,gene.end],context,true), view = span;
    const domain = [Math.min(0,span[0]),Math.max(length-1,span[1])];
    const x = value => 80 + 730 * (value-domain[0]) / Math.max(1,domain[1]-domain[0]);
    html += `<details open class="gene-track"><summary>${esc(gene.name)}${gene.name!==gene.id?' ('+esc(gene.id)+')':''} · ${esc(gene.strand)} · ${esc(gene.start)}–${esc(gene.end)} · ${esc(t('seq.overlap'))}: ${overlapBases(span,product)} bp</summary>`;
    const transcripts = gene.transcripts.length ? gene.transcripts : [{id:gene.id,segments:[]}];
    for (const transcript of transcripts) {
      let svg = `<svg class="locus-map" viewBox="0 0 900 90" role="img" aria-label="${esc(transcript.id)}"><title>${esc(transcript.id)}</title><line x1="${x(view[0])}" x2="${x(view[1])}" y1="40" y2="40" stroke="#666" stroke-width="2"/>`;
      const exons = transcript.segments.filter(segment => segment.type === 'exon');
      const overlapping = [];
      for (const [index,exon] of [...exons].sort((a,b)=>gene.strand==='-' ? b.start-a.start : a.start-b.start).entries()) {
        if (overlapBases(localSpan([exon.start,exon.end],context,true),product)) overlapping.push(index+1);
      }
      for (const segment of [...transcript.segments].sort((a,b)=>Number(a.type==='cds')-Number(b.type==='cds'))) {
        const local = localSpan([segment.start,segment.end],context,true);
        if (local[1] < local[0]) continue;
        const cds = segment.type==='cds';
        svg += `<rect x="${x(local[0])}" y="${cds?28:33}" width="${Math.max(2,x(local[1])-x(local[0]))}" height="${cds?24:14}" fill="${cds?'#846733':'#9fb9cf'}"><title>${esc(segment.type)} ${esc(segment.start)}–${esc(segment.end)}</title></rect>`;
      }
      const amplified = product;
      svg += `<rect x="${x(amplified[0])}" y="65" width="${Math.max(2,x(amplified[1])-x(amplified[0]))}" height="8" fill="#1b6255"/><text x="80" y="20">${esc(transcript.id)} · ${gene.strand === (context.strand||'+') ? '→' : '←'}</text></svg>`;
      html += svg + `<p>${esc(t('seq.exons'))}: ${esc(overlapping.join(', ') || '—')} · exon / CDS</p>`;
    }
    html += '</details>';
  }
  return html + `<p class="hint">GFF3: ${esc(annotation.source)}</p>`;
}
function sequencePage(view, start) {
  const {context,forward,reverse,pair} = view, sequence = context.sequence;
  const reverseOligo = pair && pair.reverse ? complement(pair.reverse) : '';
  start = Math.max(0,Math.min(sequence.length-1,Math.trunc(Number(start)||0)));
  const end = Math.min(sequence.length,start+1200), low=Math.min(...forward,...reverse), high=Math.max(...forward,...reverse);
  const exons = (context.annotations && context.annotations.genes || []).flatMap(gene=>gene.transcripts.flatMap(tx=>tx.segments.filter(seg=>seg.type==='exon').map(seg=>localSpan([seg.start,seg.end],context,true))));
  let html = `<p>${start+1}–${end} / ${sequence.length} bp · ${esc(t('seq.legend'))}</p><div class="sequence-table">`;
  for (let row=start;row<end;row+=60) {
    let bases=''; const last=Math.min(end,row+60);
    for (let i=row;i<last;i++) {
      const cls = i>=forward[0] && i<=forward[1] ? 'base-f' : i>=reverse[0] && i<=reverse[1] ? 'base-r' : i>=low && i<=high ? 'base-product' : '';
      const expected = pair && cls==='base-f' ? pair.forward[i-forward[0]] : cls==='base-r' ? reverseOligo[i-reverse[0]] : '';
      const changed = expected && expected.toUpperCase()!==sequence[i].toUpperCase();
      bases += `<span class="${cls}${exons.some(span=>i>=span[0]&&i<=span[1])?' base-exon':''}${changed?' base-changed':''}" title="${esc(displayCoordinate(i,context))}${changed?' '+esc(sequence[i]+' → '+expected):''}">${esc(sequence[i])}</span>`;
    }
    html += `<div class="sequence-row"><span>${esc(displayCoordinate(row,context))}</span><code>${bases}</code><span>${esc(displayCoordinate(last-1,context))}</span></div>`;
  }
  return html+'</div>';
}
function sequenceLevelView(context,forward,reverse,pair) {
  const id=sequenceViews.push({context,forward,reverse,pair})-1;
  const start=Math.max(0,Math.min(...forward,...reverse)-60);
  return `<details class="sequence-full" data-sequence-view="${id}" open><summary>${esc(t('seq.title'))}</summary><div class="sequence-nav"><label>${esc(t('seq.start'))}<input class="sequence-start" type="number" min="1" max="${context.sequence.length}" value="${start+1}"></label><button type="button" class="ghost sequence-move" data-step="0">${esc(t('seq.show'))}</button><button type="button" class="ghost sequence-move" data-step="-1200">${esc(t('seq.prev'))}</button><button type="button" class="ghost sequence-move" data-step="1200">${esc(t('seq.next'))}</button></div><div class="sequence-body">${sequencePage(sequenceViews[id],start)}</div></details>`;
}
function cleavageWindow(sequence,event,showCuts=true) {
  const start=Math.max(0,Math.min(event.site_pos-8,event.top_cut-5,event.bottom_cut-5));
  const end=Math.min(sequence.length,Math.max(event.site_pos+event.recognition.length+8,event.top_cut+5,event.bottom_cut+5));
  const lower=complement(sequence).split('').reverse().join('');
  let top='',bottom='';
  for (let i=start;i<=end;i++) {
    if (showCuts && (i===event.top_cut || i===event.bottom_cut)) {
      top += i===event.top_cut ? '<b class="cut-boundary">│</b>' : ' ';
      bottom += i===event.bottom_cut ? '<b class="cut-boundary">│</b>' : ' ';
    }
    if (i===end) break;
    const motif=i>=event.site_pos && i<event.site_pos+event.recognition.length;
    top += `<span${motif?' class="recognition-base"':''}>${esc(sequence[i])}</span>`;
    bottom += `<span${motif?' class="recognition-base"':''}>${esc(lower[i])}</span>`;
  }
  return `<p>${start+1}–${end} bp · ${esc(event.enzyme)} · ${esc(sequence.slice(event.site_pos,event.site_pos+event.recognition.length))} (${esc(event.recognition)}, ${esc(event.site_strand)}) · ${showCuts?'↑ '+event.top_cut+' / ↓ '+event.bottom_cut:esc(t('seq.compare'))}</p><pre class="cleavage-sequence">5′ ${top} 3′\n3′ ${bottom} 5′</pre>`;
}
function restrictionSequences(pair,context,derived) {
  const caps=pair.caps;
  if (derived && derived.allele_a_sequence) return derived;
  if (!derived && caps.allele_a_sequence) return caps;
  const screened=derived ? derived.specificity : pair, oligos=derived || pair;
  const forward=localSpan(screened.forward_pos,context,true), reverse=localSpan(screened.reverse_pos,context,true);
  if (!context.sequence || !forward || !reverse) return {};
  const low=Math.min(...forward,...reverse),high=Math.max(...forward,...reverse);
  let seq=context.sequence.slice(low,high+1);
  if (derived) seq=oligos.forward+seq.slice(oligos.forward.length,seq.length-oligos.reverse.length)+complement(oligos.reverse);
  const snp=caps.snp_local_index-low;
  if (!(snp>=0 && snp<seq.length) || !caps.alt_base) return {};
  return {allele_a_sequence:seq,allele_b_sequence:seq.slice(0,snp)+caps.alt_base+seq.slice(snp+1)};
}
function restrictionSitesView(result,sequences) {
  let html=typeof enzymeResultDetails === 'function' ? enzymeResultDetails(result) : '';
  html+=`<p class="hint">${esc(t('seq.cutLegend'))}</p>`;
  for (const [name,allele] of [['AA','a'],['BB','b']]) {
    const sequence=sequences[`allele_${allele}_sequence`];
    const cuts=(result[`allele_${allele}_cuts`]||[]).filter(cut=>cut.complete);
    html+=`<h5>${name}</h5>`;
    html+=!sequence ? `<p>${esc(t('map.noSequence'))}</p>` : cuts.length ? cuts.map(cut=>cleavageWindow(sequence,cut)).join('') : `<p>${esc(t('seq.noCut'))}</p>`;
    if (sequence && !cuts.length) html+=(result[`allele_${allele==='a'?'b':'a'}_cuts`]||[]).filter(cut=>cut.complete).map(cut=>cleavageWindow(sequence,cut,false)).join('');
  }
  return html;
}
function restrictionViews(pair,context) {
  const caps=pair.caps, derived=caps.dcaps && caps.dcaps.best;
  let html='';
  if (derived) html+=`<details class="restriction-detail" open><summary>dCAPS · ${esc(derived.enzyme)} · ${esc(t('seq.cutTitle'))}</summary>${restrictionSitesView(derived.digest,restrictionSequences(pair,context,derived))}</details>`;
  const cutting=(caps.natural_candidates||[]).filter(result=>[...(result.allele_a_cuts||[]),...(result.allele_b_cuts||[])].some(cut=>cut.complete));
  if (cutting.length) html+=`<details class="restriction-detail"><summary>${esc(t('seq.enzymes'))} · ${cutting.length}</summary>${cutting.map(result=>`<details><summary>${esc(result.enzyme)} · ${esc(result.recognition)}</summary>${restrictionSitesView(result,restrictionSequences(pair,context))}</details>`).join('')}</details>`;
  return html;
}
document.addEventListener('DOMContentLoaded',()=>{
  for (const mode of ['assay','sequence']) {
    const form=$(`[data-mode="${mode}"]`), label=document.createElement('label');
    label.innerHTML='<span data-i18n="seq.annotation"></span><input name="annotation_gff3" placeholder="/path/annotation.gff3.gz">';
    form.querySelector('.actions').insertAdjacentElement('beforebegin',label);
  }
  document.addEventListener('click',event=>{
    const button=event.target.closest('.sequence-move'); if (!button) return;
    const container=button.closest('[data-sequence-view]'),view=sequenceViews[Number(container.dataset.sequenceView)];
    const input=container.querySelector('.sequence-start');
    const start=Math.max(0,Math.min(view.context.sequence.length-1,Number(input.value)-1+Number(button.dataset.step)));
    input.value=start+1; container.querySelector('.sequence-body').innerHTML=sequencePage(view,start);
  });
  applyI18n();
});
