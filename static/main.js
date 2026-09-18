const $ = id => document.getElementById(id);
let last = null, inputs = null, catalog = null;
const selectedGlycans = new Set();
const modelNames = {free: 'Free neutral glycan mass', attached: 'Attached glycan residue mass'};
const fixed = (n, digits = 11) => Number(n).toFixed(digits);
const signed = n => `${n >= 0 ? '+' : '−'}${fixed(Math.abs(n))}`;
function el(tag, text, cls) {
  const n = document.createElement(tag);
  if (text !== undefined) n.textContent = text;
  if (cls) n.className = cls;
  return n;
}
function button(text, handler, cls = 'secondary') {
  const b = el('button', text, cls); b.type = 'button'; b.onclick = handler; return b;
}
function link(text, href) {
  const a = el('a', text); a.href = href; a.target = '_blank'; a.rel = 'noreferrer'; return a;
}
function status(text, error = false) { $('status').textContent = text; $('status').className = error ? 'error' : ''; }
async function api(url, options) {
  const r = await fetch(url, options), data = await r.json();
  if (!r.ok) throw Error(data.error || 'Request failed');
  return data;
}
async function analyzeBatch(submitted) {
  const lines = submitted.observations.trim().split(/\r?\n/);
  if (lines.length > 100) throw Error('Enter between 1 and 100 measurements.');
  let combined = null;
  for (let start = 0; start < lines.length; start += 10) {
    const chunk = {...submitted, observations: lines.slice(start,start+10).join('\n')};
    const result = await api('/api/analyze',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(chunk)});
    if (!combined) combined = result;
    else combined.results.push(...result.results);
    if (lines.length > 10) status(`Analyzed ${Math.min(start+10,lines.length)} of ${lines.length} measurements…`);
  }
  return combined;
}
function invalidate() {
  if (last) status('Inputs changed. Run the search again to update the candidates and calculations.');
  last = null; $('results').hidden = true; $('calculation').hidden = true;
}
function formula() {
  const attached = $('glycan_mass_type').value === 'attached';
  const positive = $('ion_mode').value === 'positive';
  $('formula').textContent = `Peptide neutral mass = Σ(residue masses) + 18.010564684 Da H₂O. Neutral glycopeptide = peptide + glycan${attached ? '' : ' − H₂O'}. Calculated m/z = neutral mass / z ${positive ? '+' : '−'} 1.007276466621. Target neutral mass = observed m/z × z ${positive ? '−' : '+'} z × 1.007276466621.`;
  $('library-convention').textContent = `Library values below are ${attached ? 'attached glycan mass increments' : 'free reducing glycan masses'} (monoisotopic). Library selections adjust automatically when you change the model.`;
  $('custom-convention').textContent = `Custom masses are interpreted as ${modelNames[$('glycan_mass_type').value].toLowerCase()} (monoisotopic). ${attached ? 'No attachment water is subtracted.' : 'One H₂O is subtracted on attachment.'} Switching mass type converts library selections, but does not change custom numbers.`;
  if (catalog) renderLibrary();
}
$('glycan_mass_type').onchange = formula;
$('ion_mode').onchange = formula;
$('form').addEventListener('input', invalidate);
$('form').addEventListener('change', invalidate);
for (const id of ['sequence', 'sites']) $(id).addEventListener('input', () => {
  if (!$('import-info').hidden) {
    $('import-info').replaceChildren(el('p', 'Sequence or sites edited after import. Check that all positions still correspond to this sequence.', 'hint'));
  }
});
formula();

$('example').onclick = async () => {
  try {
    const d = await api('/api/example');
    invalidate(); $('sequence').value = d.sequence; $('sites').value = '51,O\n432,N\n491,N\n630,N';
    $('import-info').hidden = false;
    $('import-info').replaceChildren(el('p', 'Supplied serotransferrin example: S51, N432 (complex), N491 (atypical), N630 (complex). Import P02787 for current UniProt annotations.', 'hint'));
    status('Serotransferrin loaded. Enter your measured m/z and charge.');
  } catch (e) { status(e.message, true); }
};
$('import').onclick = async () => {
  const b = $('import'); b.disabled = true; $('example').disabled = true; $('run').disabled = true;
  status('Fetching protein sequence and glycosylation annotations from UniProt…');
  try {
    const d = await api('/api/uniprot/' + encodeURIComponent($('accession').value.trim()));
    invalidate(); $('sequence').value = d.sequence; $('sites').value = d.sites_text;
    renderAnnotations(d);
    const n = d.sites.filter(s => s.type === 'N').length, o = d.sites.filter(s => s.type === 'O').length;
    status(`Imported ${d.accession}: ${n} N-linked and ${o} O-linked annotated sites.${d.sites.length ? ' Review the imported sites, then add measurements.' : ' No usable annotated sites were imported; enter known sites manually.'}`);
  } catch (e) { status(e.message, true); }
  finally { b.disabled = false; $('example').disabled = false; $('run').disabled = false; }
};
function renderAnnotations(d) {
  const box = $('import-info'); box.hidden = false; box.replaceChildren();
  box.append(link(`${d.name} · ${d.accession} ↗`, d.source));
  box.append(el('p', `UniProt annotations${d.annotation_date ? ' · updated ' + d.annotation_date : ''}. Positions refer to the imported full sequence.`, 'hint'));
  const chips = el('div', undefined, 'site-chips');
  for (const site of d.sites) chips.append(el('span', `${site.type}-linked · ${site.residue}${site.position}`, 'site-chip'));
  box.append(chips);
  if (!d.sites.length) box.append(el('p', 'No usable N/O site annotations for this sequence. This does not mean the protein is unglycosylated. Known sites can be entered manually.', 'hint'));
  const details = el('details', undefined, 'annotation-details');
  details.append(el('summary', 'Annotation descriptions & evidence'));
  for (const site of d.sites) {
    const row = el('div', undefined, 'annotation-row');
    row.append(el('strong', `${site.residue}${site.position} · ${site.type}-linked`));
    row.append(el('p', site.descriptions.join('; ')));
    if (!site.evidences.length) row.append(el('p', 'No evidence code supplied for this annotation.', 'hint'));
    for (const evidence of site.evidences) {
      const label = [evidence.evidenceCode, evidence.source, evidence.id].filter(Boolean).join(' · ');
      const p = el('p', undefined, 'hint');
      if (evidence.source === 'PubMed' && /^\d+$/.test(evidence.id || '')) p.append(link(label, `https://pubmed.ncbi.nlm.nih.gov/${evidence.id}/`));
      else p.textContent = label;
      row.append(p);
    }
    details.append(row);
  }
  if (d.sites.length) box.append(details);
  for (const item of d.skipped) box.append(el('p', `Not imported${item.position ? ' at ' + item.position : ''}: ${item.description}. ${item.reason}.`, 'hint'));
  for (const warning of d.warnings) box.append(el('p', warning, 'note'));
  box.append(el('p', 'Annotations may be experimental or inferred; sequence motifs are not treated as confirmed sites.', 'hint'));
}

function glycanMass(g) { return $('glycan_mass_type').value === 'attached' ? g.attached_mass : g.free_mass; }
function compositionText(g, full = false) {
  return Object.entries(g.composition).map(([s, count]) => `${count} × ${full ? catalog.sugars[s].name : s}`).join(' + ');
}
function svgNode(tag, attrs = {}, text) {
  const n = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [key, value] of Object.entries(attrs)) n.setAttribute(key, value);
  if (text !== undefined) n.textContent = text;
  return n;
}
function sugarShape(sugar, x, y, size = 10) {
  const data = catalog.sugars[sugar], attrs = {fill: data.color, stroke: '#3c574a', 'stroke-width': 1.1};
  let shape;
  if (data.shape === 'circle') shape = svgNode('circle', {...attrs, cx: x, cy: y, r: size});
  if (data.shape === 'square') shape = svgNode('rect', {...attrs, x: x-size, y: y-size, width: size*2, height: size*2});
  if (data.shape === 'diamond') shape = svgNode('polygon', {...attrs, points: `${x},${y-size-2} ${x+size+2},${y} ${x},${y+size+2} ${x-size-2},${y}`});
  if (data.shape === 'triangle') shape = svgNode('polygon', {...attrs, points: `${x},${y-size-2} ${x+size+2},${y+size} ${x-size-2},${y+size}`});
  shape.append(svgNode('title', {}, `${sugar}: ${data.name}`)); return shape;
}
function glycanDiagram(g, compact = false) {
  const nodes = g.nodes, children = nodes.map(() => []), depths = nodes.map(() => 0), ys = [];
  nodes.forEach((n, i) => { if (n.parent !== null) { children[n.parent].push(i); depths[i] = depths[n.parent]+1; } });
  let leaves = 0;
  function position(i) {
    if (!children[i].length) ys[i] = 34 + leaves++ * 70;
    else { children[i].forEach(position); ys[i] = children[i].reduce((v, child) => v + ys[child], 0) / children[i].length; }
  }
  position(0);
  const maxDepth = Math.max(...depths), xs = depths.map(d => 30 + (maxDepth-d)*84), width = maxDepth*84+(compact ? 60 : 170), height = Math.max(95, leaves*70);
  const svg = svgNode('svg', {viewBox: `0 0 ${width} ${height}`, role: 'img', 'aria-label': `${g.name}, representative structure. ${compositionText(g)}.`, class: compact ? 'glycan-svg compact' : 'glycan-svg'});
  svg.append(svgNode('title', {}, `${g.name} — ${compositionText(g, true)}`));
  nodes.forEach((n, i) => {
    if (n.parent === null) return;
    const p = n.parent;
    svg.append(svgNode('line', {x1: xs[i], y1: ys[i], x2: xs[p], y2: ys[p], stroke: '#879788', 'stroke-width': 1.8}));
    svg.append(svgNode('text', {x: (xs[i]+xs[p])/2, y: (ys[i]+ys[p])/2-6, 'text-anchor': 'middle', class: 'bond-label'}, n.linkage));
  });
  if (!compact) {
    svg.append(svgNode('line', {x1: xs[0]+10, y1: ys[0], x2: xs[0]+87, y2: ys[0], stroke: '#879788', 'stroke-width': 1.8}));
    const proteinLink = g.type === 'N' ? 'β–N' : g.id === 'o_glcnac' ? 'β–O' : 'α–O';
    svg.append(svgNode('text', {x: xs[0]+46, y: ys[0]-8, 'text-anchor': 'middle', class: 'bond-label'}, proteinLink));
    svg.append(svgNode('text', {x: xs[0]+93, y: ys[0]+4, class: 'protein-label'}, g.type === 'N' ? 'Asn' : 'Ser/Thr'));
  }
  nodes.forEach((n,i) => { svg.append(sugarShape(n.sugar, xs[i], ys[i])); svg.append(svgNode('text', {x: xs[i], y: ys[i]+26, 'text-anchor': 'middle', class: 'sugar-label'}, n.sugar)); });
  return svg;
}
function openStructure(g) {
  $('structure-title').textContent = g.name;
  const box = $('structure-content'); box.replaceChildren(glycanDiagram(g));
  box.append(el('p', 'Representative structure; this is not a structure assignment from your measured mass.', 'hint'));
  if (g.note) box.append(el('p', g.note, 'hint'));
  const list = el('ul', undefined, 'subunit-list');
  for (const [s,count] of Object.entries(g.composition)) list.append(el('li', `${count} × ${s} — ${catalog.sugars[s].name} (${fixed(catalog.sugars[s].residue_mass,4)} Da per sugar residue)`));
  box.append(list);
  box.append(el('p', `Attached increment: ${fixed(g.attached_mass)} Da. Free glycan: ${fixed(g.free_mass)} Da = attached increment + ${fixed(catalog.water)} Da H₂O.`, 'note'));
  box.append(el('p', 'Free masses use ExPASy GlycanMass (underivatized, monoisotopic), including its reducing end. Attached increments subtract the calculation model’s water mass.', 'hint'));
  box.append(link('ExPASy GlycanMass ↗',catalog.mass_source));
  box.append(link('Glycan family reference ↗',g.source));
  $('structure-dialog').showModal();
}
$('close-structure').onclick = () => $('structure-dialog').close();
$('structure-dialog').addEventListener('click', e => { if (e.target === $('structure-dialog')) { const r = e.target.getBoundingClientRect(); if (e.clientX < r.left || e.clientX > r.right || e.clientY < r.top || e.clientY > r.bottom) e.target.close(); } });
function updateLibrarySelection() {
  $('library-count').textContent = `${selectedGlycans.size} SELECTED`;
  $('select-all-glycans').disabled = !catalog || selectedGlycans.size === catalog.glycans.length;
  $('clear-glycans').disabled = selectedGlycans.size === 0;
}
$('select-all-glycans').onclick = () => {
  if (!catalog) return;
  catalog.glycans.forEach(g => selectedGlycans.add(g.id));
  invalidate(); renderLibrary();
};
$('clear-glycans').onclick = () => {
  selectedGlycans.clear(); invalidate(); renderLibrary();
};
function renderLibrary() {
  const box = $('glycan-library');
  const openTypes = new Set([...box.querySelectorAll('details[open]')].map(d => d.dataset.type));
  box.replaceChildren();
  for (const kind of ['O','N']) {
    const entries = catalog.glycans.filter(g => g.type === kind), details = el('details',undefined,'glycan-group');
    details.dataset.type = kind; details.open = openTypes.has(kind);
    details.append(el('summary', `Browse ${entries.length} ${kind}-linked glycans`));
    const grid = el('div',undefined,'glycan-grid');
    for (const g of entries) {
      const card = el('div',undefined,'glycan-card'), label = el('label',undefined,'glycan-select');
      const checkbox = el('input'); checkbox.type = 'checkbox'; checkbox.value = g.id; checkbox.checked = selectedGlycans.has(g.id);
      checkbox.setAttribute('aria-label', `Include ${g.name}`);
      checkbox.onchange = () => { if (checkbox.checked) selectedGlycans.add(g.id); else selectedGlycans.delete(g.id); updateLibrarySelection(); };
      label.append(checkbox,el('span',g.name)); card.append(label,glycanDiagram(g,true));
      card.append(el('p',`${fixed(glycanMass(g),7)} Da`,'glycan-mass'),el('p',compositionText(g),'hint'));
      card.append(button('Structure & subunit names ↗',()=>openStructure(g),'structure-button'));
      grid.append(card);
    }
    details.append(grid); box.append(details);
  }
  updateLibrarySelection();
}
async function loadCatalog() {
  try {
    catalog = await api('/api/glycans'); renderLibrary();
    for (const [s,data] of Object.entries(catalog.sugars)) {
      const item = el('span'), icon = svgNode('svg',{viewBox:'0 0 30 30','aria-hidden':'true'});
      icon.append(sugarShape(s,15,15,7)); item.append(icon,document.createTextNode(`${s} · ${data.name}`)); $('sugar-legend').append(item);
    }
  } catch (e) {
    $('glycan-library').replaceChildren(el('p','The library could not load. Custom masses are still available.','error'),button('Retry library',loadCatalog));
  }
}
loadCatalog();
function parseMassFile(text) {
  const tokens = text.replace(/^\uFEFF/,'').trim().split(/[,;\s]+/).filter(Boolean).map(t=>t.replace(/^"|"$/g,''));
  if (tokens.length && /^(mass|mass_da|mw)$/i.test(tokens[0])) tokens.shift();
  if (!tokens.length) throw Error('The file has no masses.');
  if (tokens.some(t => !/^(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(t) || !Number.isFinite(Number(t)) || Number(t) <= 0)) throw Error('Use positive numeric masses only, separated by commas or newlines. Optional header: mass, mass_da, or mw.');
  return tokens.map(Number);
}
for (const kind of ['O','N']) $('file_'+kind).onchange = async e => {
  const file = e.target.files[0]; if (!file) return;
  try {
    if (file.size > 50000) throw Error('Mass files must be smaller than 50 KB.');
    const imported = parseMassFile(await file.text());
    const old = $('glycans_'+kind).value.trim(), existing = old ? parseMassFile(old) : [];
    const masses = [...new Set([...existing,...imported])];
    if (masses.length > 20) throw Error('Use at most 20 custom masses per type. Clear unwanted masses before importing.');
    $('glycans_'+kind).value = masses.join(', '); invalidate();
    $('mass-import-status').textContent = `Imported ${imported.length} ${kind}-linked mass values from ${file.name}; ${masses.length} unique values ready.`;
    $('mass-import-status').className = 'hint';
  } catch (error) { $('mass-import-status').textContent = error.message; $('mass-import-status').className = 'error'; }
  e.target.value = '';
};

$('form').onsubmit = async e => {
  e.preventDefault(); $('run').disabled = true; $('import').disabled = true; invalidate(); status('Searching contiguous peptides…');
  // Snapshot inputs so later edits cannot change the meaning of returned results.
  const submitted = {...Object.fromEntries(new FormData($('form'))), selected_glycans: [...selectedGlycans]};
  try {
    const result = await analyzeBatch(submitted);
    const current = {...Object.fromEntries(new FormData($('form'))), selected_glycans: [...selectedGlycans]};
    if (JSON.stringify(current) !== JSON.stringify(submitted)) { status('Inputs changed during the search. Run again to use your latest settings.'); return; }
    inputs = submitted; last = result; render();
    status('Search complete. Closest 25 candidates shown per measurement. Select “View calculation” to inspect the mass arithmetic.');
  } catch (error) { status(error.message,true); }
  finally { $('run').disabled = false; $('import').disabled = false; }
};
function peptideCode(r) {
  const code = el('code'), offset = r.site-r.start;
  code.append(document.createTextNode(r.sequence.slice(0,offset)),el('mark',r.sequence[offset]),document.createTextNode(r.sequence.slice(offset+1)));
  return code;
}
function render() {
  const container = $('result-list'); container.replaceChildren();
  $('summary').textContent = last.candidate_count.toLocaleString()+' candidate assignments evaluated';
  last.results.forEach((result,ri) => {
    const block = el('div',undefined,'result-block'); block.append(el('h2',`m/z ${result.mz} · charge ${result.charge}${last.ion_mode === 'positive' ? '+' : '−'}`));
    block.append(el('p',`${result.match_count} within tolerance · neutral target ${fixed(result.target_mass)} Da${result.match_count===0 ? ' · No match within tolerance. These are the nearest alternatives.' : ''}`));
    const wrap = el('div',undefined,'table-wrap'),table=el('table'),head=el('tr');
    for (const title of ['#','Peptide tag / position','Site','Glycan Da','Predicted m/z','|Δm/z|','Signed ppm','Tolerance','Calculation']) head.append(el('th',title));
    const thead=el('thead'); thead.append(head); table.append(thead); const body=el('tbody');
    result.rows.forEach((r,i) => {
      const tr=el('tr'); tr.append(el('td',String(i+1)));
      const cell=el('td'); cell.append(peptideCode(r),el('p',`${r.start}–${r.end}`)); tr.append(cell);
      tr.append(el('td',`${r.type} · ${r.site}`));
      const glycan = el('td',fixed(r.glycan)); glycan.append(el('p',r.glycan_info.names.join(' / ') || 'Custom mass', 'glycan-result-name')); tr.append(glycan);
      for (const text of [fixed(r.predicted_mz),fixed(r.absolute_delta_mz),fixed(r.error_ppm,2)]) tr.append(el('td',text));
      const state=el('td'); state.append(el('span',r.within?'Within':'Outside',r.within?'yes':'no')); tr.append(state);
      const action=el('td'); action.append(button('View calculation',()=>selectCalculation(ri,i,true),'structure-button')); tr.append(action); body.append(tr);
    });
    table.append(body); wrap.append(table); block.append(wrap); container.append(block);
  });
  $('results').hidden = false;
  $('calc-measurement').replaceChildren();
  last.results.forEach((r,i) => { const option=el('option',`${i+1}. m/z ${r.mz} · ${r.charge}${last.ion_mode === 'positive' ? '+' : '−'}`); option.value=i; $('calc-measurement').append(option); });
  selectCalculation(0,0,false); $('results').scrollIntoView({behavior:'smooth',block:'start'});
}
$('calc-measurement').onchange = () => selectCalculation(Number($('calc-measurement').value),0,false);
$('calc-candidate').onchange = () => renderCalculation(Number($('calc-measurement').value),Number($('calc-candidate').value));
function selectCalculation(ri,ci,scroll) {
  if (!last) return;
  $('calc-measurement').value=ri; $('calc-candidate').replaceChildren();
  last.results[ri].rows.forEach((r,i) => { const option=el('option',`#${i+1} · ${r.start}–${r.end} · ${r.type}${r.site} · ${fixed(r.glycan,4)} Da · ${r.within?'within':'outside'}`); option.value=i; $('calc-candidate').append(option); });
  $('calc-candidate').value=ci; renderCalculation(ri,ci); $('calculation').hidden=false;
  if (scroll) $('calculation').scrollIntoView({behavior:'smooth',block:'start'});
}
function renderCalculation(ri,ci) {
  if (!last) return;
  const result=last.results[ri], r=result.rows[ci], box=$('calculation-body'); box.replaceChildren();
  if (!r) { box.append(el('p','No candidates in this peptide length range.')); return; }
  const b=r.breakdown;
  box.append(el('p',modelNames[last.glycan_mass_type] + ' · ' + last.ion_mode + ' ions','eyebrow'));
  const tag=el('div',undefined,'calculation-peptide'); tag.append(peptideCode(r)); box.append(tag);
  box.append(el('p',`${r.sequence.length} amino acids · positions ${r.start}–${r.end} · glycan attached at ${r.sequence[r.site-r.start]}${r.site} · ${r.within?'within':'outside'} tolerance`, 'hint'));
  const metrics=el('div',undefined,'mass-metrics');
  for (const [label,value] of [['Target neutral mass',`${fixed(result.target_mass)} Da`],['Calculated neutral mass',`${fixed(r.mass)} Da`],['|Δm/z|',fixed(r.absolute_delta_mz)],['Signed ppm error',`${fixed(r.error_ppm,3)} ppm`]]) {
    const metric=el('div'); metric.append(el('span',label),el('strong',value)); metrics.append(metric);
  }
  box.append(metrics,el('h3','1. Amino acids used in this peptide'));
  const wrap=el('div',undefined,'table-wrap'),table=el('table'),head=el('tr');
  const aaLabel='Residue mass';
  for (const text of ['AA','Name','Count',`${aaLabel} / AA`,'H₂O subtracted / AA','Contribution (Da)']) head.append(el('th',text));
  const thead=el('thead'); thead.append(head); table.append(thead); const tbody=el('tbody');
  for (const aa of b.amino_acids) { const row=el('tr'); for (const value of [aa.aa,aa.name,aa.count,aa.input_mass_text,aa.water_per_aa_text,`${aa.count} × ${aa.residue_mass_text} = ${aa.subtotal_text}`]) row.append(el('td',String(value))); tbody.append(row); }
  table.append(tbody); wrap.append(table); box.append(wrap);
  box.append(el('p',`Sum of residue contributions = ${fixed(b.aa_total)} Da. No water is subtracted between residues.`, 'equation'));
  box.append(el('h3','2. Add the glycan & account for water'));
  const glycanBox=el('div',undefined,'calculation-glycan');
  glycanBox.append(el('strong',r.glycan_info.names.join(' / ') || 'Custom glycan mass'));
  if (r.glycan_info.ids.length && catalog) {
    for (const id of r.glycan_info.ids) {
      const g=catalog.glycans.find(g=>g.id===id); if (!g) continue;
      const ref=el('div',undefined,'calculation-glycan-reference'); ref.append(glycanDiagram(g,true),el('p',compositionText(g,true),'hint'),button(`View ${g.name} structure`,()=>openStructure(g),'structure-button')); glycanBox.append(ref);
    }
    glycanBox.append(el('p','These are the selected reference structures, not a unique identification from mass.','hint'));
  } else glycanBox.append(el('p','No structure was supplied for this custom mass. Mass alone cannot establish its sugar subunits or connectivity.','hint'));
  box.append(glycanBox);
  box.append(el('p',`Peptide neutral mass = ${fixed(b.aa_total)} residues + ${fixed(b.terminal_water)} terminal H₂O = ${fixed(b.peptide_mass)} Da.`, 'equation'));
  box.append(el('p',`Attached glycan mass = ${fixed(r.glycan)}${b.attachment_water ? ` − ${fixed(b.attachment_water)} H₂O` : ' (already attached; no water subtraction)'} = ${fixed(b.attached_glycan_mass)} Da.`, 'equation'));
  box.append(el('p',`Neutral glycopeptide = ${fixed(b.peptide_mass)} + ${fixed(b.attached_glycan_mass)} = ${fixed(r.mass)} Da.`, 'equation'));
  box.append(el('h3','3. Compare with your measured target'));
  const positive=last.ion_mode==='positive', sign=positive?'+':'−', inverse=positive?'−':'+';
  box.append(el('p',`${positive?'Positive [M + zH]ᶻ⁺':'Negative [M − zH]ᶻ⁻'} · charge magnitude z = ${result.charge}`, 'hint'));
  box.append(el('p',`Target neutral mass = ${result.mz} × ${result.charge} ${inverse} (${result.charge} × ${last.proton_mass}) = ${fixed(result.target_mass)} Da`, 'equation'));
  box.append(el('p',`Calculated m/z = ${fixed(r.mass)} / ${result.charge} ${sign} ${last.proton_mass} = ${fixed(r.predicted_mz)}`, 'equation'));
  box.append(el('p',`Signed Δm/z = ${fixed(r.predicted_mz)} − ${result.mz} = ${signed(r.signed_delta_mz)}. Absolute Δm/z = ${fixed(r.absolute_delta_mz)}.`, 'equation'));
  box.append(el('p',`Signed ppm error = (${signed(r.signed_delta_mz)} / ${result.mz}) × 10⁶ = ${signed(r.error_ppm)} ppm. ΔMass (neutral) = ${signed(r.error_da)} Da.`, 'equation'));
  box.append(el('p',`Tolerance: ${inputs.tolerance} ${inputs.unit === 'Da' ? 'Da on neutral mass' : inputs.unit === 'mz' ? 'm/z' : 'ppm'}. Positive signed error means calculated m/z is above observed m/z. Calculations use unrounded values.`, 'hint'));

}
$('download').onclick = () => {
  if (!last) return;
  const keys=['sequence','start','end','site','type','glycan','neutral_mass','predicted_mz','absolute_delta_mz','signed_delta_mz','error_da','error_ppm','absolute_ppm_error','within'];
  const totals=['aa_input_total','aa_water_loss','aa_total','terminal_water','attachment_water','peptide_mass','attached_glycan_mass'];
  const config=['glycan_mass_type','ion_mode','tolerance','unit','min_length','max_length','sites','glycans_O','glycans_N'];
  const rows=[['observed_mz','charge','target_neutral_mass',...keys,'glycan_names','glycan_ids',...totals,...config,'selected_glycans','protein_sequence']];
  for (const result of last.results) for (const r of result.rows) rows.push([result.mz,result.charge,result.target_mass,...keys.map(k=>r[k]),r.glycan_info.names.join('; '),r.glycan_info.ids.join('; '),...totals.map(k=>r.breakdown[k]),...config.map(k=>inputs[k]),last.selected_glycans.join('; '),last.sequence]);
  const quote=x=>'"'+String(x).replaceAll('"','""')+'"';
  const blob=new Blob([rows.map(r=>r.map(quote).join(',')).join('\r\n')],{type:'text/csv;charset=utf-8'}),url=URL.createObjectURL(blob),a=el('a');
  a.href=url; a.download='glycopeptide-candidates.csv'; a.click(); setTimeout(()=>URL.revokeObjectURL(url),1000);
};
