// Review imported measurements before changing the search inputs.
(() => {
  let peaks = [], selected = new Set(), page = 0, importVersion = 0;
  const pageSize = 50;
  const message = (text, error = false) => {
    $('peak-status').textContent = text;
    $('peak-status').className = error ? 'error' : 'hint';
  };
  const eligible = peak => peak.charge !== null &&
    (!peak.polarity || peak.polarity === $('peak-ion-mode').value);
  const onPage = () => peaks.slice(page * pageSize, (page + 1) * pageSize);
  function updateSelection() {
    $('peak-add').textContent = `Add ${selected.size} selected to measurements`;
    $('peak-add').disabled = selected.size === 0 || !$('peak-ion-mode').value;
    $('peak-clear').disabled = selected.size === 0;
    $('peak-select').disabled = !$('peak-ion-mode').value || !onPage().some(eligible);
  }
  function renderPeaks() {
    const body = $('peak-rows'); body.replaceChildren();
    onPage().forEach((peak, offset) => {
      const index = page * pageSize + offset;
      const row = el('tr'), cell = el('td'), check = el('input');
      check.type = 'checkbox'; check.checked = selected.has(index);
      check.disabled = !eligible(peak) || !$('peak-ion-mode').value;
      check.setAttribute('aria-label', `Select peak ${index + 1}, m/z ${peak.mz}`);
      check.onchange = () => {
        if (check.checked && selected.size >= 100) {
          check.checked = false;
          message('Select at most 100 peaks per search.', true); return;
        }
        if (check.checked) selected.add(index); else selected.delete(index);
        updateSelection();
      };
      cell.append(check); row.append(cell);
      const z = peak.charge === null ? 'Unknown' : `${peak.charge}${peak.polarity === 'positive' ? '+' : peak.polarity === 'negative' ? '−' : ''}`;
      for (const value of [peak.mz, z, peak.intensity ?? '—', peak.scan ?? '—']) row.append(el('td', String(value)));
      body.append(row);
    });
    $('peak-page').textContent = `${page * pageSize + 1}–${Math.min((page + 1) * pageSize, peaks.length)} of ${peaks.length} peaks`;
    $('peak-prev').disabled = page === 0;
    $('peak-next').disabled = (page + 1) * pageSize >= peaks.length;
    updateSelection();
  }
  $('raw-file').onchange = () => {
    const file = $('raw-file').files[0];
    $('raw-status').textContent = file
      ? (!/\.raw$/i.test(file.name) ? 'Choose a Thermo .raw file.' : `${file.name} selected. Analysis is unavailable until a RAW backend is connected. This file has not been uploaded or analyzed. You can import a peak table from Qual alongside it.`)
      : 'Selecting a file does not upload or analyze it. Connect a RAW analysis backend before importing binary .raw data.';
  };
  $('peak-file').onchange = async () => {
    const version = ++importVersion, file = $('peak-file').files[0];
    peaks = []; selected.clear(); page = 0;
    $('peak-review').hidden = true; $('peak-ion-mode').value = '';
    if (!file) { message(''); return; }
    if (!/\.(csv|tsv|txt)$/i.test(file.name)) { message('Choose a CSV, TSV, or TXT peak table exported from Qual.', true); return; }
    if (!file.size || file.size > 2000000) { message('Peak tables must be between 1 byte and 2 MB.', true); return; }
    message(`Reading ${file.name}…`);
    try {
      const data = await api('/api/import-peaks', {
        method: 'POST', headers: {'Content-Type': 'application/octet-stream'}, body: file
      });
      if (version !== importVersion) return;
      peaks = data.peaks;
      if (data.polarities.length === 1) $('peak-ion-mode').value = data.polarities[0];
      $('peak-provenance').textContent = `${file.name} · ${data.source}. ${data.polarities.length > 1 ? 'Mixed polarities: select one ion mode at a time.' : 'Confirm ion mode before adding measurements.'}`;
      $('peak-review').hidden = false; renderPeaks();
      message(`Loaded ${peaks.length} peaks. ${data.unknown_charge_count} have unknown charge and cannot be selected. Your existing measurements have not changed.`);
    } catch (error) { if (version === importVersion) message(error.message, true); }
  };
  $('peak-ion-mode').onchange = () => { selected.clear(); renderPeaks(); };
  $('peak-prev').onclick = () => { if (page > 0) { page--; renderPeaks(); } };
  $('peak-next').onclick = () => { if ((page + 1) * pageSize < peaks.length) { page++; renderPeaks(); } };
  $('peak-clear').onclick = () => { selected.clear(); renderPeaks(); };
  $('peak-select').onclick = () => {
    onPage().forEach((peak, offset) => {
      if (eligible(peak) && selected.size < 100) selected.add(page * pageSize + offset);
    });
    renderPeaks();
    if (selected.size === 100) message('100 peaks selected — the maximum per search.');
  };
  $('peak-add').onclick = () => {
    const mode = $('peak-ion-mode').value;
    if (!selected.size || !mode) return;
    const existing = $('observations').value.trim().split(/\r?\n/).filter(line => line.trim());
    if (existing.length && $('ion_mode').value !== mode) {
      message('The imported ion mode differs from your existing measurements. Use a separate search for each polarity; clear the measurements before adding this table.', true); return;
    }
    const seen = new Set(existing.map(line => {
      const values = line.trim().split(/[,\s]+/);
      return values.length === 2 ? `${Number(values[0])},${Number(values[1])}` : line;
    }));
    const added = [];
    for (const index of [...selected].sort((a, b) => a - b)) {
      const peak = peaks[index];
      if (!eligible(peak)) continue;
      const line = `${peak.mz},${peak.charge}`;
      if (!seen.has(line)) { seen.add(line); added.push(line); }
    }
    if (existing.length + added.length > 100) { message('Adding these peaks would exceed 100 measurements. Select fewer peaks or clear existing measurements.', true); return; }
    if (!added.length) { message('Those m/z–charge pairs are already in your measurements.'); return; }
    $('observations').value = [...existing, ...added].join('\n');
    $('ion_mode').value = mode; formula(); invalidate();
    message(`Added ${added.length} unique m/z–charge pairs in ${mode} ion mode. Review the measurements, then run the candidate search.`);
    $('observations').focus();
  };
})();
