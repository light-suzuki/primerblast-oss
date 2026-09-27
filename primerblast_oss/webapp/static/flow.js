/* Keep the shared search controls beside the input that they apply to. */
Object.assign(I18N.ja, {
  'flow.input': '1. 増幅したい配列・領域を入力',
  'flow.checkInput': '1. 調べたいプライマー配列を入力',
  'flow.dbInput': 'ゲノムFASTAと保存先を入力',
  'flow.db': '2. オフターゲットを検索するゲノムを選択',
  'flow.run': '3. 入力した条件で計算する',
  'flow.type': '入力欄', 'flow.select': '選択欄', 'flow.file': 'ファイルを選択',
  'flow.demo': '練習用の配列を入力する（デモ）',
  'flow.demoHint': 'デモは入力内容を置き換えます。自分のデータを使う場合は、下の欄に入力してください。',
  'flow.demoLoaded': 'デモ配列を入力しました。検索先を選んでから計算してください。',
  'flow.options': '詳細設定を変更する（通常は初期値のまま）',
  'flow.help': '必要なデータ・用語・よくある質問',
  'flow.purpose': '目的を選び直す',
  'flow.dbHint': '複数選択できます。色の付いたチェック済みのゲノムが検索対象です。',
  'flow.placeholder': '薄い例の文字は入力済みの値ではありません。'
});
Object.assign(I18N.en, {
  'flow.input': '1. Enter the sequence or region to amplify',
  'flow.checkInput': '1. Enter primers to check',
  'flow.dbInput': 'Enter a genome FASTA and output path',
  'flow.db': '2. Select genomes to search for off-targets',
  'flow.run': '3. Calculate with these inputs',
  'flow.type': 'Type here', 'flow.select': 'Choose an option', 'flow.file': 'Choose a file',
  'flow.demo': 'Fill with a practice sequence (demo)',
  'flow.demoHint': 'The demo replaces your inputs. To use your own data, fill in the fields below.',
  'flow.demoLoaded': 'Demo sequence loaded. Select search genomes before calculating.',
  'flow.options': 'Change advanced settings (defaults are usually sufficient)',
  'flow.help': 'Required data, terminology and common questions',
  'flow.purpose': 'Choose a different task',
  'flow.dbHint': 'Multiple selections are allowed. Checked, highlighted genomes will be searched.',
  'flow.placeholder': 'Faint example text is a hint, not an entered value.'
});

document.addEventListener('DOMContentLoaded', () => {
  const db = $('.db-panel');
  const dbTitle = document.createElement('h3'); dbTitle.dataset.i18n = 'flow.db'; db.prepend(dbTitle);
  const dbHint = document.createElement('p'); dbHint.className = 'flow-hint'; dbHint.dataset.i18n = 'flow.dbHint'; dbTitle.after(dbHint);
  const guide = $('.guide');
  const grid = $('.mode-grid', guide);
  const help = document.createElement('details'); help.className = 'guide-help';
  const summary = document.createElement('summary'); summary.dataset.i18n = 'flow.help'; help.append(summary);
  let next = grid.nextElementSibling;
  while (next && !next.classList.contains('guide-cta')) { const following = next.nextElementSibling; help.append(next); next = following; }
  grid.after(help);
  const store = $('.project-store'); if (store) store.open = false;

  $$('form[data-mode]').forEach(form => {
    const first = $('fieldset', form); if (!first) return;
    const heading = document.createElement('h3'); heading.className = 'flow-heading'; heading.dataset.i18n = form.dataset.mode === 'check' ? 'flow.checkInput' : form.dataset.mode === 'makedb' ? 'flow.dbInput' : 'flow.input'; first.before(heading);
    const hint = document.createElement('p'); hint.className = 'flow-hint'; hint.dataset.i18n = 'flow.placeholder'; first.before(hint);
    const example = $('.example-btn', form);
    if (example) {
      const demo = document.createElement('details'); demo.className = 'demo-choice';
      const caption = document.createElement('summary'); caption.dataset.i18n = 'flow.demo';
      const explanation = document.createElement('p'); explanation.dataset.i18n = 'flow.demoHint';
      example.dataset.i18n = 'flow.demo'; demo.append(caption, explanation, example); first.before(demo);
      const message = document.createElement('p'); message.className = 'flow-demo-message'; message.hidden = true; message.setAttribute('role', 'status'); message.dataset.i18n = 'flow.demoLoaded'; demo.after(message);
      example.addEventListener('click', () => { message.hidden = false; demo.open = false; });
      form.addEventListener('reset', () => { message.hidden = true; });
    }
    const actions = $('.actions', form);
    const options = document.createElement('details'); options.className = 'flow-options adv';
    const caption = document.createElement('summary'); caption.dataset.i18n = 'flow.options'; options.append(caption);
    let sibling = first.nextElementSibling;
    while (sibling && sibling !== actions) { const following = sibling.nextElementSibling; options.append(sibling); sibling = following; }
    if (options.children.length > 1) actions.before(options);
    const run = $('.run', form); if (run && form.dataset.mode !== 'makedb') run.dataset.i18n = 'flow.run';
  });

  $$('label').forEach(label => {
    const control = $('input:not([type="hidden"]),textarea,select', label);
    if (!control || ['checkbox', 'radio'].includes(control.type)) return;
    const tag = document.createElement('small'); tag.className = 'control-kind';
    tag.dataset.i18n = control.tagName === 'SELECT' ? 'flow.select' : control.type === 'file' ? 'flow.file' : 'flow.type';
    label.prepend(tag);
  });
  const originalShowTab = showTab;
  showTab = name => {
    originalShowTab(name);
    const form = $(`form[data-mode="${name}"]`);
    if (form && name !== 'makedb') $('.actions', form).before(db);
    else $('#tabs').before(db);
  };
  $('.tab-guide [data-i18n]').dataset.i18n = 'flow.purpose';
  const syncChips = () => {
    $$('button', db).forEach(button => { button.type = 'button'; });
    $$('.db-chip', db).forEach(chip => {
      const selected = chip.classList.contains('on');
      chip.setAttribute('aria-pressed', String(selected));
    });
  };
  new MutationObserver(syncChips).observe($('#db-list'), { childList: true, subtree: true, attributes: true, attributeFilter: ['class'] });
  syncChips(); applyI18n();
});
