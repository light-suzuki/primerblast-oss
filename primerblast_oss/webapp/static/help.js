/* Purpose-first documentation. Shortcuts only select inputs; they never run jobs. */
const experimentHelp = {
  ja: {
    tab: '説明', title: '実験の目的から、使う機能を選ぶ',
    intro: 'リファレンスを選ぶ → 対象を指定する → 検索DBを確認する → 実行する。下のボタンは入力画面を開くだけです。',
    recipes: '実験ユースケース', functions: '各タブでできること', input: '用意するもの', output: '得られるもの', open: '入力画面へ',
    read: ['結果の読み方', 'F/Rはプライマーの結合位置と向き、bpは増幅産物の長さです。遺伝子注釈があればエクソン・イントロンと増幅範囲を確認できます。オフターゲットは各産物のbpと検索状態を確認してください。探索が打ち切られた結果は「他の産物がない」ことの証明にはなりません。'],
    data: ['データの保持と共有', 'プロジェクト保存はこのブラウザ内の保存です。バックアップや別の環境への移動にはJSONを書き出し、読み込みで復元します。保存・書き出しには配列やローカルパスが含まれ得ます。共有前に内容を確認してください。ゲノムの登録と解析はローカルサーバーで行います。将来のエージェント連携でも、参照の登録は外部送信の許可を意味しません。'],
    wet: ['計算結果を実験へ持ち込むとき', '仮想ゲルは予測図です。PCRの増幅、酵素消化、ヘテロ接合体の判別は実験で確認する必要があります。小さい断片や近いサイズのバンドが実際に見分けられるか、親系統・陰性対照と併せて確認してください。'],
    recipesList: [
      ['遺伝子をPCRで増幅したい', 'リファレンスと遺伝子IDを選び、遺伝子全体・CDSなど対象を指定します。', 'design', 'gene'],
      ['遺伝子の配列を確かめたい', '遺伝子IDから重複するPCR断片を設計し、カバーする範囲を確認します。', 'sequence', 'gene'],
      ['特定の領域を読みたい', '染色体と座標で領域を指定し、シーケンス用の重複PCRを設計します。', 'sequence', 'interval'],
      ['SNPで遺伝子型を判定したい', 'SNP位置と代替塩基からCAPS・dCAPS・AS-PCR候補と予測バンドを比較します。', 'assay', 'snp'],
      ['手持ちのプライマーを確認したい', 'F/R配列を入力し、選んだゲノムDBで増幅産物とオフターゲットを調べます。', 'check', ''],
      ['QTL領域にマーカーを配置したい', '領域と配置数または間隔から、複数位置のPCR候補を設計します。多型の確認は別途必要です。', 'markers', '']
    ],
    functionsList: [
      ['PCR設計', '遺伝子ID＋ゲノム＋注釈、またはDNA配列。特異性確認には検索DB。', 'プライマー候補、増幅サイズ、結合位置、特異性の結果。', 'design'],
      ['PCR確認', '既存のF/Rプライマー配列と検索DB。', '予測増幅産物の位置・bpとオフターゲット。', 'check'],
      ['シーケンス', '遺伝子ID、ゲノム座標、またはDNA配列。', '重複PCRの候補とカバー範囲。実測の読長は保証しません。', 'sequence'],
      ['ゲル・CAPS', '遺伝子、領域、またはSNPと代替塩基。必要に応じてVCF。', '制限酵素の認識配列・切断位置、断片長、仮想ゲル、利用可能なマーカー候補。', 'assay'],
      ['QTLマーカー', 'ゲノム座標とマーカー数または間隔。', '領域内に配置するPCR候補。QTL原因遺伝子の推定機能ではありません。', 'markers'],
      ['ゲノムの準備', 'ローカルのゲノムFASTA。', '特異性検索に使うBLAST DB。遺伝子IDの検索には別途注釈が必要です。', 'makedb'],
      ['タイルPCR（プロ向け）', 'DNA配列と断片サイズ・重複長。', '長い配列を複数の重複PCRに分割する候補。通常のシーケンス設計はシーケンスタブから始められます。', 'tile']
    ]
  },
  en: {
    tab: 'Help', title: 'Choose a workflow for your experiment',
    intro: 'Choose a reference → specify a target → review search databases → run. The buttons below only open input screens.',
    recipes: 'Experimental use cases', functions: 'What each tab does', input: 'Required inputs', output: 'Outputs', open: 'Open inputs',
    read: ['Reading results', 'F/R indicate primer binding positions and directions; bp is the product length. With gene annotation, inspect exons, introns and the amplified region. Review each off-target product size and the search status. A truncated search does not establish the absence of other products.'],
    data: ['Saving and sharing data', 'Saved projects stay in this browser. Export JSON for backup or transfer, and import it to restore a project. Saved and exported data may contain sequences and local paths: review before sharing. Reference registration and analysis use the local server. Future agent integration must not treat reference registration as permission to transmit data externally.'],
    wet: ['Taking predictions to the bench', 'The virtual gel is a prediction. Validate amplification, digestion and heterozygote discrimination experimentally. Check whether small fragments or closely spaced bands can be resolved, alongside parental samples and negative controls.'],
    recipesList: [
      ['Amplify a gene', 'Choose a reference and gene ID, then specify the whole gene, CDS or another feature.', 'design', 'gene'],
      ['Check a gene sequence', 'Design overlapping PCR products from a gene ID and review their coverage.', 'sequence', 'gene'],
      ['Sequence a genomic region', 'Specify a chromosome interval and design overlapping PCR products for sequencing.', 'sequence', 'interval'],
      ['Genotype a SNP', 'Provide a SNP position and alternate base to compare CAPS, dCAPS and AS-PCR candidates and predicted bands.', 'assay', 'snp'],
      ['Check existing primers', 'Enter F/R sequences to inspect products and off-targets against selected genome databases.', 'check', ''],
      ['Place markers across a QTL interval', 'Choose an interval and marker count or spacing to design PCR candidates. Polymorphisms require separate confirmation.', 'markers', '']
    ],
    functionsList: [
      ['PCR design', 'Gene ID with genome and annotation, or DNA sequence. Search databases for specificity checks.', 'Primer candidates, product sizes, binding positions and specificity results.', 'design'],
      ['PCR check', 'Existing F/R sequences and search databases.', 'Predicted product positions, sizes and off-targets.', 'check'],
      ['Sequencing', 'Gene ID, genomic interval or DNA sequence.', 'Overlapping PCR candidates and coverage. Actual sequencing read length is not guaranteed.', 'sequence'],
      ['Gel / CAPS', 'Gene, interval, or SNP with alternate base; optional VCF.', 'Enzyme recognition sequences, cleavage positions, fragment sizes, virtual gels and available marker candidates.', 'assay'],
      ['QTL markers', 'Genomic interval and marker count or spacing.', 'PCR candidates distributed across the interval; this does not infer causal QTL genes.', 'markers'],
      ['Prepare genome', 'Local genome FASTA.', 'BLAST database for specificity searches. Gene-ID lookup separately requires annotation.', 'makedb'],
      ['Tile PCR (advanced)', 'DNA sequence, product size and overlap.', 'Candidates dividing a long sequence into overlapping PCR products. For routine sequencing, start with the Sequencing tab.', 'tile']
    ]
  }
};
Object.entries(experimentHelp).forEach(([lang, words]) => {
  ['tab', 'title', 'intro', 'recipes', 'functions'].forEach(key => { I18N[lang][`help.${key}`] = words[key]; });
  ['read', 'data', 'wet'].forEach(key => {
    I18N[lang][`help.${key}.title`] = words[key][0];
    I18N[lang][`help.${key}.body`] = words[key][1];
  });
});
function renderExperimentHelp() {
  const words = experimentHelp[LANG] || experimentHelp.en;
  $('#experiment-recipes').innerHTML = words.recipesList.map(([title, body, mode, source]) =>
    `<article class="help-card"><h4>${esc(title)}</h4><p>${esc(body)}</p><button type="button" data-help-mode="${mode}" data-help-source="${source}">${esc(words.open)}</button></article>`).join('');
  $('#function-help').innerHTML = words.functionsList.map(([title, input, output, mode]) =>
    `<article class="help-card"><h4>${esc(title)}</h4><dl><dt>${esc(words.input)}</dt><dd>${esc(input)}</dd><dt>${esc(words.output)}</dt><dd>${esc(output)}</dd></dl><button type="button" data-help-mode="${mode}">${esc(words.open)}</button></article>`).join('');
}
document.addEventListener('DOMContentLoaded', () => {
  const previousApply = applyI18n;
  applyI18n = () => { previousApply(); renderExperimentHelp(); };
  const previousShow = showTab;
  showTab = name => {
    previousShow(name);
    const help = name === 'help';
    if (help) $('.db-panel').hidden = true;
    // Preserve results and saved projects while displaying documentation.
    $('#results').classList.toggle('help-hidden', help);
    const store = $('.project-store'); if (store) store.classList.toggle('help-hidden', help);
  };
  $('[data-panel="help"]').addEventListener('click', event => {
    const button = event.target.closest('[data-help-mode]'); if (!button) return;
    const mode = button.dataset.helpMode; showTab(mode);
    const source = button.dataset.helpSource;
    const selector = mode === 'design' ? '#design-src' : mode === 'sequence' ? '#sequence-src' : mode === 'assay' ? '#assay-kind' : null;
    const control = selector && $(selector);
    if (source && control) { control.value = source; control.dispatchEvent(new Event('change')); }
  });
  applyI18n();
});
