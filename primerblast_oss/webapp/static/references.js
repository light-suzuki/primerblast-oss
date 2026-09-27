Object.assign(I18N.ja, {
  'ref.title': '先にリファレンスゲノムを選ぶ', 'ref.manual': '手動でファイルパスを指定する',
  'ref.hint': '選ぶとFASTA・GFF3・検索DBが入力されます。入力した遺伝子IDは保持します。',
  'ref.applied': '選んだゲノムのファイルパスを入力しました。検索DBも選択済みです。',
  'ref.noAnnotation': 'このゲノムには対応するGFF3が登録されていません。遺伝子IDで設計するには注釈を指定してください。',
  'ref.noDb': '検索DBは未登録です。下で検索先を選択してください。',
  'ref.failed': 'ゲノム一覧を取得できませんでした。手動入力で利用できます。',
  'ref.pcrPurpose': '1. 何をPCRで増幅したい？', 'ref.pcrGene': '遺伝子IDから設計',
  'ref.pcrSequence': '手元の配列から設計',
  'ref.idBuilder': '遺伝子IDの数字部分を入力', 'ref.chromosome': '染色体番号',
  'ref.geneNumber': '遺伝子番号', 'ref.builderHint': '固定の文字は入力不要です。番号の先頭の0は自動で補います。名前・別名を使う場合は下のID欄に直接入力できます。',
  'ref.geneHint': '全角でも入力できます。注釈に登録されたID・名前・別名に対応します。'
});
Object.assign(I18N.en, {
  'ref.title': 'Choose a reference genome first', 'ref.manual': 'Enter file paths manually',
  'ref.hint': 'Fills FASTA, GFF3 and search database paths. Your gene ID is retained.',
  'ref.applied': 'Reference file paths filled. The search database is selected too.',
  'ref.noAnnotation': 'No matching GFF3 is registered. Supply an annotation to design by gene ID.',
  'ref.noDb': 'No search database is registered. Select one below.',
  'ref.failed': 'Could not load reference genomes. Manual input remains available.',
  'ref.pcrPurpose': '1. What do you want to amplify by PCR?', 'ref.pcrGene': 'Design from a gene ID',
  'ref.pcrSequence': 'Design from your own sequence',
  'ref.idBuilder': 'Enter the numeric parts of the gene ID', 'ref.chromosome': 'Chromosome number',
  'ref.geneNumber': 'Gene number', 'ref.builderHint': 'Fixed text is supplied. Leading zeros are added automatically. Enter a name or alias directly in the ID field below.',
  'ref.geneHint': 'Full-width characters are accepted. Use an ID, name or alias registered in the annotation.'
});

function composeReferenceGeneId(format, chromosome, number) {
  const chrom = String(chromosome).normalize('NFKC').trim();
  const digits = String(number).normalize('NFKC').trim();
  if (!format || !Number.isInteger(format.digits) || format.digits < 1 || format.digits > 12 ||
      !/^\d+$/.test(chrom) || !/^\d+$/.test(digits) || digits.length > format.digits ||
      !Array.isArray(format.chromosomes) || !format.chromosomes.includes(chrom)) return '';
  return `${format.prefix}${chrom}${format.separator}${digits.padStart(format.digits, '0')}`;
}

function syncGeneBuilder(builder, format) {
  const valid = format && typeof format.prefix === 'string' && typeof format.separator === 'string' &&
    Number.isInteger(format.digits) && format.digits >= 1 && format.digits <= 12 &&
    Array.isArray(format.chromosomes) && format.chromosomes.length && format.chromosomes.every(chrom => typeof chrom === 'string' && /^\d+$/.test(chrom));
  builder.format = valid ? format : null; builder.box.hidden = !valid;
  if (!valid) return;
  $('.id-prefix', builder.box).textContent = format.prefix;
  $('.id-separator', builder.box).textContent = format.separator;
  builder.number.maxLength = format.digits; builder.number.placeholder = '0'.repeat(format.digits);
  builder.chromosome.replaceChildren(...format.chromosomes.map(chrom => { const option = document.createElement('option'); option.value = chrom; option.textContent = chrom; return option; }));
  const raw = builder.gene.value.normalize('NFKC').trim().replace(/^gene\s*:\s*/i, '');
  builder.number.value = '';
  for (const chrom of format.chromosomes) {
    const head = `${format.prefix}${chrom}${format.separator}`;
    if (raw.toLowerCase().startsWith(head.toLowerCase())) {
      const number = raw.slice(head.length);
      if (/^\d+$/.test(number) && number.length === format.digits) {
        builder.chromosome.value = chrom; builder.number.value = number; break;
      }
    }
  }
}

document.addEventListener('DOMContentLoaded', async () => {
  const designSource = $('#design-src');
  const designForm = designSource.form;
  const purpose = document.createElement('fieldset'); purpose.className = 'experiment-choices pcr-choices';
  purpose.innerHTML = '<legend data-i18n="ref.pcrPurpose"></legend>';
  const choices = [];
  [['gene', 'ref.pcrGene'], ['sequence', 'ref.pcrSequence']].forEach(([value, key]) => {
    const label = document.createElement('label'); label.className = 'experiment-choice';
    const input = document.createElement('input'); input.type = 'radio'; input.name = '_pcrPurpose'; input.dataset.uiControl = 'true'; input.value = value;
    const caption = document.createElement('span'); caption.dataset.i18n = key;
    label.append(input, caption); purpose.append(label); choices.push(input);
    input.addEventListener('change', () => { designSource.value = value; designSource.dispatchEvent(new Event('change')); });
  });
  $('.mode-intro', designForm).after(purpose); $('.flow-heading', designForm).hidden = true;
  designSource.closest('label').hidden = true;
  Array.from(designSource.options).forEach(option => { option.defaultSelected = option.value === 'gene'; });
  const syncPurpose = () => choices.forEach(input => { input.checked = input.value === designSource.value; input.closest('label').classList.toggle('chosen', input.checked); });
  designSource.addEventListener('change', syncPurpose);
  designForm.addEventListener('reset', () => setTimeout(syncPurpose, 0));
  designSource.value = 'gene'; designSource.dispatchEvent(new Event('change')); syncPurpose();
  const feature = $('select[name="gene_feature"]', designForm);
  const featureLabel = $('span', feature.closest('label'));
  featureLabel.dataset.i18n = 'ref.pcrFeature';
  I18N.ja['ref.pcrFeature'] = '遺伝子のどの範囲を増幅する？';
  I18N.en['ref.pcrFeature'] = 'Which gene span do you want to amplify?';
  Array.from(feature.options).forEach(option => { const value = option.value; option.value = value; option.dataset.i18n = `flow.feature.${value}`; });

  const selectors = [];
  $$('form[data-mode]').forEach(form => {
    if (!$('input[name="genome"]', form)) return;
    const block = document.createElement('section'); block.className = 'reference-picker';
    block.innerHTML = '<label><small class="control-kind" data-i18n="flow.select"></small><span data-i18n="ref.title"></span><select><option value="" data-i18n="ref.manual"></option></select></label><p class="hint" data-i18n="ref.hint"></p><p role="status" class="reference-status"></p>';
    $('.mode-intro', form).after(block);
    const select = $('select', block); select.id = `reference-${form.dataset.mode}`;
    const gene = $('input[name="gene"]', form);
    let builder = null;
    if (gene) {
      const box = document.createElement('section'); box.className = 'gene-id-builder'; box.hidden = true;
      const group = gene.closest('[data-src],[data-seq-src],[data-kind]');
      if (group) { Object.assign(box.dataset, group.dataset); box.style.display = group.style.display; }
      box.innerHTML = '<strong data-i18n="ref.idBuilder"></strong><div class="gene-id-parts"><code class="id-prefix"></code><label><small data-i18n="ref.chromosome"></small><select></select></label><code class="id-separator"></code><label><small data-i18n="ref.geneNumber"></small><input inputmode="numeric" autocomplete="off"></label></div><p class="hint" data-i18n="ref.builderHint"></p>';
      gene.closest('label').before(box);
      const chromosome = $('select', box), number = $('input', box);
      builder = {box, chromosome, number, gene, format: null};
      const write = () => {
        number.value = number.value.normalize('NFKC');
        gene.value = composeReferenceGeneId(builder.format, chromosome.value, number.value);
      };
      chromosome.addEventListener('change', write); number.addEventListener('input', write);
      gene.addEventListener('input', () => syncGeneBuilder(builder, builder.format));
    }
    selectors.push({form, select, status: $('.reference-status', block), builder});
    $$('input[name="gene"]', form).forEach(input => {
      const hint = document.createElement('em'); hint.className = 'hint'; hint.dataset.i18n = 'ref.geneHint'; input.after(hint);
    });
  });
  applyI18n();
  let profiles = [];
  let autoDatabase = '';
  const syncReferences = () => selectors.forEach(({form, select, builder}) => {
    const genome = form.elements.namedItem('genome').value;
    const annotation = form.elements.namedItem('gff3');
    const index = profiles.findIndex(profile => profile.available && profile.genome === genome && (!annotation || profile.gff3 === annotation.value));
    select.value = index < 0 ? '' : String(index);
    if (builder) syncGeneBuilder(builder, index < 0 ? null : profiles[index].gene_id_format);
  });
  const previousShowTab = showTab;
  showTab = name => { previousShowTab(name); syncReferences(); };
  selectors.forEach(({form, select, status, builder}) => {
    select.addEventListener('change', () => {
      if (select.value === '') { status.textContent = ''; if (builder) syncGeneBuilder(builder, null); return; }
      const profile = profiles[Number(select.value)]; if (!profile || !profile.available) return;
      if (builder) syncGeneBuilder(builder, profile.gene_id_format);
      ['genome', 'gff3', 'annotation_gff3'].forEach(name => {
        const input = form.elements.namedItem(name);
        if (input) input.value = name === 'genome' ? profile.genome : profile.gff3;
      });
      const keys = [profile.database ? 'ref.applied' : 'ref.noDb'];
      if (!profile.gff3) keys.push('ref.noAnnotation');
      status.replaceChildren(...keys.map(key => { const line = document.createElement('span'); line.dataset.i18n = key; return line; }));
      selectedDbs = selectedDbs.filter(db => db !== autoDatabase);
      if (profile.database) {
        selectedDbs = [profile.database, ...selectedDbs.filter(db => db !== profile.database)];
        customDbs = [...new Set([...customDbs, profile.database])];
      }
      const mapping = $('#db-genomes');
      const lines = mapping.value.split(/\r?\n/).filter(line => line.trim() && selectedDbs.includes(line.split('=')[0].trim()) && line.split('=')[0].trim() !== profile.database);
      if (profile.database) {
        lines.push(`${profile.database}=${profile.genome}`);
      }
      mapping.value = lines.join('\n'); autoDatabase = profile.database;
      loadDatabases();
      applyI18n();
    });
    $$('input[name="genome"],input[name="gff3"]', form).forEach(input => input.addEventListener('change', () => { syncReferences(); status.textContent = ''; }));
    form.addEventListener('reset', () => setTimeout(() => { syncReferences(); status.textContent = ''; }, 0));
  });
  try {
    const response = await fetch('/api/references'); if (!response.ok) throw new Error();
    const catalog = await response.json(); profiles = catalog.references || [];
    selectors.forEach(({select, status}) => {
      profiles.forEach((profile, index) => {
        const option = document.createElement('option'); option.value = index;
        option.textContent = profile.name + (profile.available ? '' : ` (${profile.missing.join(', ')})`);
        option.disabled = !profile.available; select.append(option);
      });
      if (catalog.warnings && catalog.warnings.length) status.textContent = catalog.warnings.join('\n');
    });
    syncReferences();
  } catch (_) {
    selectors.forEach(({status}) => { status.dataset.i18n = 'ref.failed'; }); applyI18n();
  }
});
