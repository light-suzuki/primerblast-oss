/* Local, searchable enzyme reference. Browsing never launches an analysis. */
const enzymeWords = {
  ja: {
    tab:'制限酵素', title:'制限酵素の名前・認識配列・切断パターン',
    intro:'同じ配列を認識する酵素も、切断位置やメチル化感受性・反応条件が違うことがあります。同じ切断の酵素と、切断が異なる酵素を分けて表示します。',
    search:'酵素名・確認済み製品名・認識配列で検索（全角も可）',
    legend:'「|」が切断位置。上は5′→3′、下は3′→5′です。外側のNは隣接する塩基、認識配列内のNは任意の塩基を示します。座標は認識配列の先頭からの境界位置（0始まり）。',
    prev:'前の30件', next:'次の30件', products:'確認済み製品名', same:'同じ認識配列・同じ切断', different:'同じ認識配列・異なる切断',
    unknown:'切断位置が不明', unknownRelated:'同じ認識配列・切断位置不明', supported:'単一切断モデルで配列予測可能', unsupported:'通常PCRの自動候補から除外（切断不明・複数切断・特殊基質）',
    commercial:'入手可能の記録あり（収録時点）', unavailable:'供給元の記録なし（収録時点）',
    conditions:'同じ切断でも反応条件は同一とは限りません。メチル化、複数認識部位の必要性、温度、バッファー、末端付近の切断効率は使用する製品の資料で確認してください。製品名の収録も全メーカー・全製品の網羅を保証しません。',
    source:'収録元', count:'件', matched:'検索結果', error:'酵素一覧を読み込めませんでした。ローカルサーバーへの接続を確認してください。',
    blunt:'平滑末端', over5:'5′突出末端', over3:'3′突出末端', offset:'上鎖 / 下鎖の切断境界'
  },
  en: {
    tab:'Enzymes', title:'Restriction enzyme names, recognition sites and cleavage patterns',
    intro:'Enzymes recognizing the same sequence may differ in cleavage, methylation sensitivity and reaction requirements. Same-cut and different-cut enzymes are listed separately.',
    search:'Search enzyme name, verified product name or recognition sequence (full-width accepted)',
    legend:'The | marks cleavage. Top: 5′→3′; bottom: 3′→5′. Flanking N represents adjacent bases; N within the recognition site means any base. Offsets are zero-based boundaries from the recognition sequence start.',
    prev:'Previous 30', next:'Next 30', products:'Verified product names', same:'Same recognition and cleavage', different:'Same recognition, different cleavage',
    unknown:'Cleavage positions unknown', unknownRelated:'Same recognition, unknown cleavage', supported:'Sequence prediction with a single-cut model', unsupported:'Excluded from ordinary PCR automatic candidates (unknown/multiple cuts or special substrates)',
    commercial:'Supplier recorded in snapshot', unavailable:'No supplier recorded in snapshot',
    conditions:'Identical cleavage does not imply identical reaction conditions. Check product documentation for methylation, multiple-site requirements, temperature, buffer and end-proximity effects. Product names do not cover every vendor or product.',
    source:'Source', count:'records', matched:'Matches', error:'Could not load enzymes. Check the local server connection.',
    blunt:'Blunt end', over5:'5′ overhang', over3:'3′ overhang', offset:'Top / bottom cleavage boundaries'
  }
};
Object.entries(enzymeWords).forEach(([lang, words]) => Object.entries(words).forEach(([key,value]) => { I18N[lang][`enz.${key}`]=value; }));
let enzymeCatalog = null, enzymePage = 0;
function enzymePatternHTML(info) {
  if (!info || !info.pattern) return `<p>${esc(t('enz.unknown'))}</p>`;
  return `<pre class="cleavage-sequence">${esc(info.pattern.top)}\n${esc(info.pattern.bottom)}</pre>` +
    (info.cuts || []).map(([top,bottom]) => `<p>${esc(t('enz.offset'))}: ${top} / ${bottom} · ${esc(t(top === bottom ? 'enz.blunt' : top < bottom ? 'enz.over5' : 'enz.over3'))}${top === bottom ? '' : ` (${Math.abs(top-bottom)} nt)`}</p>`).join('');
}
function enzymeRelationshipsHTML(info) {
  return [['products',info.product_names],['same',info.same_cut_enzymes],['different',info.different_cut_enzymes],['unknownRelated',info.unknown_cut_enzymes]].filter(([,names]) => names && names.length).map(([key,names]) => `<p><strong>${esc(t('enz.'+key))}</strong>: ${names.map(esc).join(', ')}</p>`).join('');
}
function enzymeResultDetails(result) {
  const info = result.pattern ? {...result,cuts:[[result.top_cut_offset,result.bottom_cut_offset]]} : enzymeCatalog && enzymeCatalog.enzymes.find(row => row.name === result.enzyme && row.recognition === result.recognition && row.cuts.length === 1 && row.cuts[0][0] === result.top_cut_offset && row.cuts[0][1] === result.bottom_cut_offset);
  return info ? enzymePatternHTML(info) + enzymeRelationshipsHTML(info) : '';
}
function renderEnzymeCatalog() {
  if (!enzymeCatalog) return;
  const query = $('#enzyme-search').value.normalize('NFKC').replace(/\s+/g,'').toLowerCase();
  const rows = enzymeCatalog.enzymes.filter(row => [row.name,row.recognition,...(row.product_names||[]),...(row.same_cut_enzymes||[]),...(row.different_cut_enzymes||[])].some(name => name.toLowerCase().includes(query)));
  const start = enzymePage * 30;
  $('#enzyme-status').textContent = `${t('enz.source')}: ${enzymeCatalog.release} / REBASE ${enzymeCatalog.rebase_version} · ${enzymeCatalog.enzymes.length} ${t('enz.count')} · ${t('enz.matched')}: ${rows.length} (${rows.length ? start+1 : 0}–${Math.min(start+30,rows.length)})`;
  $('#enzyme-records').innerHTML = rows.slice(start,start+30).map(row => `<article class="help-card"><h4>${esc(row.name)} · ${esc(row.recognition)}</h4>${enzymePatternHTML(row)}${enzymeRelationshipsHTML(row)}<p>${esc(t(row.prediction_supported ? 'enz.supported' : 'enz.unsupported'))}</p><p>${esc(t(row.commercial_in_snapshot ? 'enz.commercial' : 'enz.unavailable'))}</p><a href="${esc(row.reference_url)}" target="_blank" rel="noopener noreferrer">REBASE · ${esc(row.name)}</a></article>`).join('');
  $('#enzyme-prev').disabled = enzymePage === 0;
  $('#enzyme-next').disabled = start+30 >= rows.length;
}
document.addEventListener('DOMContentLoaded', () => {
  $('#enzyme-search').addEventListener('input', () => { enzymePage=0; renderEnzymeCatalog(); });
  $('#enzyme-prev').addEventListener('click', () => { enzymePage--; renderEnzymeCatalog(); });
  $('#enzyme-next').addEventListener('click', () => { enzymePage++; renderEnzymeCatalog(); });
  const previousApply = applyI18n;
  applyI18n = () => { previousApply(); renderEnzymeCatalog(); };
  const previousShow = showTab;
  showTab = name => {
    previousShow(name);
    if (name === 'enzymes') {
      $('.db-panel').hidden = true; $('#results').classList.add('help-hidden');
      const store=$('.project-store'); if (store) store.classList.add('help-hidden');
    }
  };
  fetch('/api/enzymes').then(response => { if (!response.ok) throw new Error(); return response.json(); }).then(data => {
    enzymeCatalog=data; renderEnzymeCatalog();
  }).catch(() => { $('#enzyme-status').textContent=t('enz.error'); });
  applyI18n();
});
