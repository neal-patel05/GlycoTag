const $ = id => document.getElementById(id);
let last = null, inputs = null, catalog = null;
let lastGlycanAdducts = null;
const selectedGlycans = new Set();
const selectedGlycosmos = new Map(), selectedAdducts = new Set();
let allGlycosmos = false;
const excludedGlycosmos = new Set();
const glycosmosRecords = new Map();
const linkageLabels = {all:'all glycans', N:'N-linked glycans', O:'O-linked glycans', both:'N- and O-linked glycans'};
let adductCatalog = null, glycosmosPageData = null, glycosmosRequest = 0;
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
    else {
      combined.results.push(...result.results);
      combined.glycosmos_image_ids = [...new Set([...(combined.glycosmos_image_ids || []), ...(result.glycosmos_image_ids || [])])];
    }
    if (lines.length > 10) status(`Analyzed ${Math.min(start+10,lines.length)} of ${lines.length} measurements…`);
  }
  return combined;
}
function invalidate() {
  if (last) status('Inputs changed. Run the search again to update the candidates and calculations.');
  last = null; $('results').hidden = true; $('calculation').hidden = true;
  lastGlycanAdducts = null; $('glycan-adduct-results').hidden = true;
}
function formula() {
  const attached = $('glycan_mass_type').value === 'attached';
  const positive = $('ion_mode').value === 'positive';
  $('formula').textContent = `Peptide neutral mass = Σ(residue masses) + 18.010564684 Da H₂O. Neutral glycopeptide = peptide + glycan${attached ? '' : ' − H₂O'}. Calculated m/z = neutral mass / z ${positive ? '+' : '−'} 1.007276466621. Target neutral mass = observed m/z × z ${positive ? '−' : '+'} z × 1.007276466621.`;
  if ($('adduct-mode').value === 'selected') $('formula').textContent = `Peptide neutral mass = Σ(residue masses) + H₂O. Neutral glycopeptide = peptide + glycan${attached ? '' : ' − H₂O'}. Each selected adduct uses m/z = (n × neutral mass + shift) / z; target neutral mass = (observed m/z × z − shift) / n. Only ${positive ? 'positive' : 'negative'} adducts matching the measurement charge are searched.`;
  $('library-convention').textContent = `Library values below are ${attached ? 'attached glycan mass increments' : 'free reducing glycan masses'} (monoisotopic). Library selections adjust automatically when you change the model.`;
  $('custom-convention').textContent = `Custom masses are interpreted as ${modelNames[$('glycan_mass_type').value].toLowerCase()} (monoisotopic). ${attached ? 'No attachment water is subtracted.' : 'One H₂O is subtracted on attachment.'} Switching mass type converts library selections, but does not change custom numbers.`;
  if (catalog) renderLibrary();
  if (adductCatalog) renderAdducts();
  if (glycosmosPageData) renderGlycosmos();
  renderGlycosmosSelection();
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
  const submitted = searchInputs();
  try {
    const result = await analyzeBatch(submitted);
    const current = searchInputs();
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
  if (last.glycosmos_coverage) {
    const c = last.glycosmos_coverage;
    container.append(el('p', `GlyCosmos: ${c.selected.toLocaleString()} selected; ${c.supported.toLocaleString()} with supported masses; ${c.skipped.toLocaleString()} skipped because no supported ExPASy mass is available.`, 'note'));
    container.append(el('p', `${c.matched.toLocaleString()} GlyCosmos records eligible for these protein sites. ${c.linkage_skipped.toLocaleString()} additional records skipped: ${c.unclassified.toLocaleString()} without N/O linkage annotations and ${c.incompatible.toLocaleString()} without a compatible enabled site type. N-linked references are matched only to N sites; O-linked references only to S/T sites.`, 'note'));
    if (c.selected) container.append(el('p','Structures show the selected GlyCosmos references. Different structures can share the same mass; a mass match does not uniquely identify a structure.','hint'));
  }
  if (last.skipped_sites?.length) container.append(el('p',`Sites skipped because no compatible glycan masses were selected: ${last.skipped_sites.map(s=>`${s.type} · ${s.site}`).join(', ')}.`, 'note'));
  last.results.forEach((result,ri) => {
    const block = el('div',undefined,'result-block'); block.append(el('h2',`m/z ${result.mz} · charge ${result.charge}${last.ion_mode === 'positive' ? '+' : '−'}`));
    block.append(el('p',`${result.match_count} within tolerance · ${result.target_mass === null ? 'neutral targets depend on the adduct' : 'neutral target ' + fixed(result.target_mass) + ' Da'}${result.match_count===0 && result.rows.length ? ' · No match within tolerance. These are the nearest alternatives.' : ''}`));
    if (result.note) block.append(el('p', result.note, 'note'));
    const wrap = el('div',undefined,'table-wrap'),table=el('table'),head=el('tr');
    for (const title of ['#','Peptide tag / position','Site','Glycan Da','Adduct','Predicted m/z','|Δm/z|','Signed ppm','Tolerance','Calculation']) head.append(el('th',title));
    const thead=el('thead'); thead.append(head); table.append(thead); const body=el('tbody');
    result.rows.forEach((r,i) => {
      const tr=el('tr'); tr.append(el('td',String(i+1)));
      const cell=el('td'); cell.append(peptideCode(r),el('p',`${r.start}–${r.end}`)); tr.append(cell);
      tr.append(el('td',`${r.type} · ${r.site}`));
      const glycan = el('td',fixed(r.glycan)); glycan.append(el('p',r.glycan_info.names.join(' / ') || 'Custom mass', 'glycan-result-name'));
      const sourceIds = r.glycan_info.ids.filter(id=>/^G\d{5}[A-Z]{2}$/.test(id));
      if (sourceIds.length) {
        glycan.append(glycosmosReference(sourceIds[0]));
        if (sourceIds.length > 1) {
          const more=el('details',undefined,'result-structures');
          more.append(el('summary',`Show ${sourceIds.length-1} more reference structures`));
          for (const id of sourceIds.slice(1)) more.append(glycosmosReference(id));
          glycan.append(more);
        }
      }
      if (r.glycan_info.library_count > r.glycan_info.ids.length) glycan.append(el('p', `${r.glycan_info.library_count.toLocaleString()} selected references share this mass; first 20 shown.`, 'hint'));
      tr.append(glycan);
      tr.append(el('td',r.adduct));
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
  last.results[ri].rows.forEach((r,i) => { const option=el('option',`#${i+1} · ${r.start}–${r.end} · ${r.type}${r.site} · ${fixed(r.glycan,4)} Da · ${r.adduct} · ${r.within?'within':'outside'}`); option.value=i; $('calc-candidate').append(option); });
  $('calc-candidate').value=ci; renderCalculation(ri,ci); $('calculation').hidden=false;
  if (scroll) $('calculation').scrollIntoView({behavior:'smooth',block:'start'});
}
function renderCalculation(ri,ci) {
  if (!last) return;
  const result=last.results[ri], r=result.rows[ci], box=$('calculation-body'); box.replaceChildren();
  if (!r) { box.append(el('p',result.note || 'No candidates in this peptide length range.')); return; }
  const b=r.breakdown;
  box.append(el('p',modelNames[last.glycan_mass_type] + ' · ' + last.ion_mode + ' ions','eyebrow'));
  const tag=el('div',undefined,'calculation-peptide'); tag.append(peptideCode(r)); box.append(tag);
  box.append(el('p',`${r.sequence.length} amino acids · positions ${r.start}–${r.end} · glycan attached at ${r.sequence[r.site-r.start]}${r.site} · ${r.within?'within':'outside'} tolerance`, 'hint'));
  const metrics=el('div',undefined,'mass-metrics');
  for (const [label,value] of [['Target neutral mass',`${fixed(r.target_neutral_mass)} Da`],['Calculated neutral mass',`${fixed(r.mass)} Da`],['|Δm/z|',fixed(r.absolute_delta_mz)],['Signed ppm error',`${fixed(r.error_ppm,3)} ppm`]]) {
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
      const g=catalog.glycans.find(g=>g.id===id);
      if (!g) { if (/^G\d{5}[A-Z]{2}$/.test(id)) glycanBox.append(glycosmosReference(id)); continue; }
      const ref=el('div',undefined,'calculation-glycan-reference'); ref.append(glycanDiagram(g,true),el('p',compositionText(g,true),'hint'),button(`View ${g.name} structure`,()=>openStructure(g),'structure-button')); glycanBox.append(ref);
    }
    glycanBox.append(el('p','These are the selected reference structures, not a unique identification from mass.','hint'));
  } else glycanBox.append(el('p','No structure was supplied for this custom mass. Mass alone cannot establish its sugar subunits or connectivity.','hint'));
  box.append(glycanBox);
  box.append(el('p',`Peptide neutral mass = ${fixed(b.aa_total)} residues + ${fixed(b.terminal_water)} terminal H₂O = ${fixed(b.peptide_mass)} Da.`, 'equation'));
  box.append(el('p',`Attached glycan mass = ${fixed(r.glycan)}${b.attachment_water ? ` − ${fixed(b.attachment_water)} H₂O` : ' (already attached; no water subtraction)'} = ${fixed(b.attached_glycan_mass)} Da.`, 'equation'));
  box.append(el('p',`Neutral glycopeptide = ${fixed(b.peptide_mass)} + ${fixed(b.attached_glycan_mass)} = ${fixed(r.mass)} Da.`, 'equation'));
  box.append(el('h3','3. Compare with your measured target'));
  box.append(el('p',`${last.ion_mode==='positive'?'Positive':'Negative'} ${r.adduct} · charge magnitude z = ${result.charge} · multimer count n = ${r.multimer}`, 'hint'));
  box.append(el('p',`Target neutral mass = (${result.mz} × ${result.charge} − (${signed(r.adduct_mass_shift)})) / ${r.multimer} = ${fixed(r.target_neutral_mass)} Da`, 'equation'));
  box.append(el('p',`Calculated m/z = (${r.multimer} × ${fixed(r.mass)} + (${signed(r.adduct_mass_shift)})) / ${result.charge} = ${fixed(r.predicted_mz)}`, 'equation'));
  box.append(el('p',`Signed Δm/z = ${fixed(r.predicted_mz)} − ${result.mz} = ${signed(r.signed_delta_mz)}. Absolute Δm/z = ${fixed(r.absolute_delta_mz)}.`, 'equation'));
  box.append(el('p',`Signed ppm error = (${signed(r.signed_delta_mz)} / ${result.mz}) × 10⁶ = ${signed(r.error_ppm)} ppm. ΔMass (neutral) = ${signed(r.error_da)} Da.`, 'equation'));
  box.append(el('p',`Tolerance: ${inputs.tolerance} ${inputs.unit === 'Da' ? 'Da on neutral mass' : inputs.unit === 'mz' ? 'm/z' : 'ppm'}. Positive signed error means calculated m/z is above observed m/z. Calculations use unrounded values.`, 'hint'));

}
$('download').onclick = () => {
  if (!last) return;
  const keys=['sequence','start','end','site','type','glycan','neutral_mass','adduct','adduct_id','adduct_mass_shift','multimer','predicted_mz','absolute_delta_mz','signed_delta_mz','error_da','error_ppm','absolute_ppm_error','within'];
  const totals=['aa_input_total','aa_water_loss','aa_total','terminal_water','attachment_water','peptide_mass','attached_glycan_mass'];
  const config=['glycan_mass_type','ion_mode','tolerance','unit','min_length','max_length','sites','glycans_O','glycans_N'];
  const rows=[['observed_mz','charge','target_neutral_mass',...keys,'glycan_names','glycan_ids','glycan_reference_count',...totals,...config,'selected_glycans','glycosmos_selection','glycosmos_coverage','selected_adducts','protein_sequence']];
  for (const result of last.results) for (const r of result.rows) rows.push([result.mz,result.charge,r.target_neutral_mass,...keys.map(k=>r[k]),r.glycan_info.names.join('; '),r.glycan_info.ids.join('; '),r.glycan_info.library_count ?? r.glycan_info.ids.length,...totals.map(k=>r.breakdown[k]),...config.map(k=>inputs[k]),last.selected_glycans.join('; '),JSON.stringify(last.glycosmos_selection),JSON.stringify(last.glycosmos_coverage),last.selected_adducts?.join('; ') ?? 'proton-only',last.sequence]);
  const quote=x=>'"'+String(x).replaceAll('"','""')+'"';
  const blob=new Blob([rows.map(r=>r.map(quote).join(',')).join('\r\n')],{type:'text/csv;charset=utf-8'}),url=URL.createObjectURL(blob),a=el('a');
  a.href=url; a.download='glypeptid-candidates.csv'; a.click(); setTimeout(()=>URL.revokeObjectURL(url),1000);
};

function searchInputs() {
  const data = {...Object.fromEntries(new FormData($('form'))), selected_glycans: [...selectedGlycans],
    glycosmos_selection: {all:allGlycosmos, ids:allGlycosmos ? [] : [...selectedGlycosmos.keys()], excluded:[...excludedGlycosmos], linkage:$('glycosmos-linkage').value, species:$('glycosmos-species').value},
    glycosmos_apply_to:$('glycosmos-apply-to').value};
  if ($('adduct-mode').value === 'selected') data.selected_adducts = [...selectedAdducts];
  return data;
}
function renderAdducts() {
  const focused = document.activeElement?.closest('.adduct-choice')?.querySelector('input')?.value;
  const choosing = $('adduct-mode').value === 'selected', mode = $('ion_mode').value;
  const count = adductCatalog.adducts.filter(a => a.mode === mode && selectedAdducts.has(a.id)).length;
  $('adduct-count').textContent = choosing ? `${count} ${mode.toUpperCase()} SELECTED` : 'PROTON ONLY';
  const box = $('adduct-library'); box.replaceChildren();
  for (const polarity of ['positive', 'negative']) {
    const entries = adductCatalog.adducts.filter(a => a.mode === polarity);
    const group = el('details', undefined, 'glycan-group adduct-group'); group.dataset.mode = polarity; group.open = true;
    group.append(el('summary', `${polarity === 'positive' ? '1. Positive' : '2. Negative'} ion mode · ${entries.length} adducts${mode === polarity ? ' · active' : ''}`));
    const actions = el('div', undefined, 'library-actions');
    const all = button(`Select all ${polarity}`, () => { $('adduct-mode').value='selected'; entries.forEach(a => selectedAdducts.add(a.id)); invalidate(); formula(); });
    all.disabled = choosing && entries.every(a => selectedAdducts.has(a.id));
    const clear = button(`Clear ${polarity}`, () => { entries.forEach(a => selectedAdducts.delete(a.id)); invalidate(); renderAdducts(); });
    clear.disabled = entries.every(a => !selectedAdducts.has(a.id));
    actions.append(all, clear); group.append(actions);
    const grid = el('div', undefined, 'adduct-grid');
    for (const a of entries) {
      const label = el('label', undefined, 'glycan-select adduct-choice'), checkbox = el('input');
      checkbox.type = 'checkbox'; checkbox.value = a.id; checkbox.checked = selectedAdducts.has(a.id);
      checkbox.onchange = () => { $('adduct-mode').value='selected'; checkbox.checked ? selectedAdducts.add(a.id) : selectedAdducts.delete(a.id); invalidate(); formula(); };
      label.append(checkbox, el('span', a.name)); grid.append(label);
    }
    group.append(grid); box.append(group);
  }
  if (focused) box.querySelector(`input[value="${focused}"]`)?.focus();
}
$('adduct-mode').onchange = formula;
async function loadAdducts() {
  try {
    adductCatalog = await api('/api/adducts');
    // Familiar defaults; selecting all remains an explicit user choice.
    adductCatalog.adducts.filter(a => ['[M+H]+','[M+2H]2+','[M+3H]3+','[M-H]-','[M-2H]2-','[M-3H]3-'].includes(a.name)).forEach(a => selectedAdducts.add(a.id));
    $('adduct-note').textContent = adductCatalog.note; renderAdducts();
  } catch (e) {
    $('adduct-library').replaceChildren(el('p', 'Adduct library could not load. ' + e.message, 'error'), button('Retry adduct library', loadAdducts));
  }
}
loadAdducts();

let glycosmosSpeciesTimer;
async function loadGlycosmos(page = 1) {
  clearTimeout(glycosmosSpeciesTimer);
  const request = ++glycosmosRequest;
  const query = new URLSearchParams({q:$('glycosmos-query').value, page:String(page), supported:$('glycosmos-supported').checked ? '1' : '0', linkage:$('glycosmos-linkage').value, species:$('glycosmos-species').value});
  $('glycosmos-status').textContent = 'Searching local GlyCosmos library…';
  $('glycosmos-status').className = 'hint';
  $('glycosmos-prev').disabled = $('glycosmos-next').disabled = true;
  $('glycosmos-select-page').disabled = $('glycosmos-all').disabled = true;
  $('glycosmos-library').replaceChildren();
  try {
    const data = await api('/api/glycosmos?' + query);
    if (request !== glycosmosRequest) return;
    glycosmosPageData = data; renderGlycosmos();
  } catch (e) {
    if (request !== glycosmosRequest) return;
    glycosmosPageData = null; $('glycosmos-library').replaceChildren(); $('glycosmos-page').textContent = '';
    $('glycosmos-status').textContent = e.message; $('glycosmos-status').className = 'error';
  }
}
function glycosmosChosen(id) { return allGlycosmos ? !excludedGlycosmos.has(id) : selectedGlycosmos.has(id); }
function chooseGlycosmos(g, checked) {
  if (allGlycosmos) checked ? excludedGlycosmos.delete(g.id) : excludedGlycosmos.add(g.id);
  else checked ? selectedGlycosmos.set(g.id,g) : selectedGlycosmos.delete(g.id);
  invalidate(); renderGlycosmosSelection();
}
function structureImage(g, large = false) {
  const box=el('div',undefined,'structure-image-box'), img=el('img');
  if (!g.has_image) {
    box.classList.add('image-unavailable');
    box.append(el('p','Structure image not supplied by GlyCosmos','hint'),el('p','View source notation and subunit details below.','hint'));
    return box;
  }
  img.src=g.image_url; img.alt=`GlyCosmos structure of ${g.id}`;
  img.loading=large?'eager':'lazy'; img.decoding='async'; img.className=large?'glycosmos-image large':'glycosmos-image';
  img.referrerPolicy='no-referrer';
  if (large) img.onload=()=>{img.style.maxWidth=Math.min(img.naturalWidth*2,800)+'px';img.style.margin='auto';};
  img.onerror=()=>{
    box.replaceChildren(el('p','The structure image could not load. Refresh the page or view the GlyCosmos source record.','hint'));
  };
  box.append(img); return box;
}
function glycosmosReference(id, result = last) {
  const ref=el('div',undefined,'glycosmos-result-reference');
  ref.append(link(id,`https://glycosmos.org/glycans/${id}`),structureImage({
    id, has_image:result.glycosmos_image_ids?.includes(id), image_url:`/api/glycosmos/image/${id}`
  }));
  const status=el('p','','hint');
  const show=button('Structure & subunit names ↗',async()=>{
    show.disabled=true;status.textContent='Loading structure details…';
    try {
      let g=glycosmosRecords.get(id);
      if (!g) { g=await api(`/api/glycosmos/record/${id}`);glycosmosRecords.set(id,g); }
      openGlycosmos(g);status.textContent='';
    } catch(e) { status.textContent=e.message; }
    finally { show.disabled=false; }
  },'structure-button');
  ref.append(show,status);
  return ref;
}
function subunitList(g, expanded = false) {
  const list=el('ul',undefined,'glycosmos-subunits');
  for (const unit of g.subunits) {
    const item=el('li',`${unit.count} × ${unit.symbol} — ${unit.name}`);
    if (expanded) {
      if (unit.source_id) item.append(document.createTextNode(' · '),link('Residue reference ↗',`https://glycosmos.org/glycans/${unit.source_id}`));
      if (unit.symbol==='Unresolved residue') item.append(el('code',unit.wurcs,'iupac-text'));
    }
    list.append(item);
  }
  if (!g.subunits.length) list.append(el('li','Subunit names are not defined by this source record.'));
  return list;
}
function openGlycosmos(g) {
  $('structure-title').textContent=g.id + ' · GlyCosmos structure';
  const box=$('structure-content');box.replaceChildren(structureImage(g,true));
  box.append(el('h3',g.name),el('p',g.name_kind,'hint'));
  box.append(el('p',g.linkage_types.length ? `${g.linkage_types.map(k=>k+'-linked').join(' / ')} · GlyCosmos motif annotation` : 'No N/O motif annotation in GlyCosmos', 'hint'));
  box.append(el('h3','Reported organisms'),el('p',organismText(g,true),'glycosmos-organisms hint'));
  if (g.motifs) box.append(el('p','Source motif annotations: '+g.motifs,'hint'));
  box.append(el('h3','Subunit names'),subunitList(g,true));
  if (g.subunit_note) box.append(el('p',g.subunit_note,'note'));
  box.append(el('p',g.mass_note,'hint'));
  if(g.free_mass!==null) box.append(el('p',`Free glycan: ${fixed(g.free_mass,7)} Da · attached increment: ${fixed(g.attached_mass,7)} Da`,'note'));
  const source=el('details');source.append(el('summary','Source WURCS notation'),el('p',g.wurcs,'iupac-text'));box.append(source);
  box.append(link('GlyCosmos record ↗',g.source));
  $('structure-dialog').showModal();
}
function organismText(g, expanded = false) {
  const organisms=g.organisms || [];
  if (!organisms.length) return 'No organism annotation supplied by GlyCosmos.';
  return organisms.slice(0,expanded?organisms.length:3).map(t=>`${t.name} (${t.id})`).join('; ')
    + (!expanded && organisms.length>3 ? `; +${organisms.length-3} more in structure details` : '');
}
function renderGlycosmos() {
  const d=glycosmosPageData, box=$('glycosmos-library');box.replaceChildren();
  const options=$('glycosmos-species-options');options.replaceChildren();
  for (const taxon of d.species_options) {
    const option=el('option');option.value=taxon.name;
    option.label=`${taxon.aliases.join(' / ')}${taxon.aliases.length?' · ':''}${taxon.id} · ${taxon.count.toLocaleString()} glycans`;
    options.append(option);
  }
  const matches=d.species_matches;
  $('glycosmos-species-status').textContent=!d.species
    ? `All organisms, including unannotated records. ${d.metadata.taxonomy_annotated.toLocaleString()} glycans have source organism annotations.`
    : !matches.length ? 'No organism annotations match this name or taxonomy ID. Try a scientific name or clear the field.'
    : `Reported in: ${matches.slice(0,5).map(t=>`${t.name} (${t.id})`).join('; ')}${matches.length>5?`; and ${matches.length-5} more matching organisms. Choose a suggestion for a specific organism.`:'.'}`;
  $('glycosmos-status').textContent=`${d.total.toLocaleString()} matching records · ${d.metadata.total.toLocaleString()} total · ${d.image_count.toLocaleString()} source structure images · ${d.metadata.supported.toLocaleString()} supported ExPASy masses. GlyCosmos · CC BY 4.0 · ${d.metadata.retrieved}.`;
  const grid=el('div',undefined,'glycan-grid');
  for (const g of d.glycans) {
    glycosmosRecords.set(g.id,g);
    const card=el('article',undefined,'glycosmos-card');
    const label=el('label',undefined,'glycan-select'),checkbox=el('input');
    checkbox.type='checkbox';checkbox.value=g.id;checkbox.checked=glycosmosChosen(g.id);
    checkbox.setAttribute('aria-label',`Include ${g.id}`);
    checkbox.onchange=()=>chooseGlycosmos(g,checkbox.checked);
    label.append(checkbox,el('span',g.id));card.append(label,el('p',g.name,'glycosmos-name'),el('p',g.name_kind,'hint'));
    card.append(el('p',g.linkage_types.length ? `${g.linkage_types.map(k=>k+'-linked').join(' / ')} · source motif` : 'No N/O motif annotation', 'hint'));
    card.append(el('p',organismText(g),'glycosmos-organisms hint'));
    card.append(structureImage(g));
    card.append(subunitList(g));
    if(g.subunit_note) card.append(el('p',g.subunit_note,'hint'));
    if(g.free_mass!==null) card.append(el('p',`${fixed(glycanMass(g),7)} Da`,'glycan-mass'));
    else card.append(el('p','No supported ExPASy mass · retained in library selections','hint'));
    card.append(button('Structure & subunit names ↗',()=>openGlycosmos(g),'structure-button'));
    grid.append(card);
  }
  box.append(grid);
  $('glycosmos-page').textContent=d.pages?`Page ${d.page} of ${d.pages.toLocaleString()}`:'No matching records';
  $('glycosmos-prev').disabled=d.page<=1; $('glycosmos-next').disabled=d.page>=d.pages;
  $('glycosmos-select-page').disabled=!d.glycans.length;
  renderGlycosmosSelection();
}
function renderGlycosmosSelection() {
  const linkage=$('glycosmos-linkage').value;
  const scope=linkageLabels[linkage]+($('glycosmos-species').value.trim() ? ` matching organism “${$('glycosmos-species').value.trim()}”` : ' across all organisms');
  const total=glycosmosPageData?.linkage_counts[linkage] ?? 0;
  const count=allGlycosmos ? total-excludedGlycosmos.size : selectedGlycosmos.size;
  $('glycosmos-count').textContent=`${count.toLocaleString()} SELECTED`;
  $('glycosmos-clear').disabled=!count;
  $('glycosmos-all').disabled=!total || (allGlycosmos&&!excludedGlycosmos.size);
  $('glycosmos-all').textContent=`Select all ${total.toLocaleString()} ${linkage==='all'?'glycans':linkageLabels[linkage]}`;
  $('glycosmos-selected').replaceChildren(el('p',allGlycosmos
    ? `All ${total.toLocaleString()} records in ${scope} selected across all pages${excludedGlycosmos.size?`, except ${excludedGlycosmos.size.toLocaleString()} unchecked records`:''}. Peptide matching uses only linkage-compatible records with supported masses; skipped records are reported.`
    : `${count.toLocaleString()} records selected across pages in ${scope}. Peptide matching requires a matching N/O linkage annotation.`, 'hint'));
}
function resetGlycosmosScope() {
  clearTimeout(glycosmosSpeciesTimer);++glycosmosRequest;glycosmosPageData=null;
  allGlycosmos=false;selectedGlycosmos.clear();excludedGlycosmos.clear();
  invalidate();renderGlycosmosSelection();
  $('glycosmos-library').replaceChildren();$('glycosmos-page').textContent='';
  $('glycosmos-prev').disabled=$('glycosmos-next').disabled=$('glycosmos-select-page').disabled=true;
  $('glycosmos-status').textContent='Searching local GlyCosmos library…';
  $('glycosmos-species-status').textContent='Updating organism filter…';
}
$('glycosmos-linkage').onchange=()=>{resetGlycosmosScope();loadGlycosmos();};
$('glycosmos-species').oninput=()=>{
  resetGlycosmosScope();glycosmosSpeciesTimer=setTimeout(()=>loadGlycosmos(),300);
};
$('glycosmos-species').onkeydown=e=>{if(e.key==='Enter'){e.preventDefault();loadGlycosmos();}};
$('glycosmos-search').onclick=()=>loadGlycosmos();
$('glycosmos-query').onkeydown=e=>{if(e.key==='Enter'){e.preventDefault();loadGlycosmos();}};
$('glycosmos-supported').onchange=()=>loadGlycosmos();
$('glycosmos-prev').onclick=()=>loadGlycosmos(glycosmosPageData.page-1);
$('glycosmos-next').onclick=()=>loadGlycosmos(glycosmosPageData.page+1);
$('glycosmos-all').onclick=()=>{allGlycosmos=true;selectedGlycosmos.clear();excludedGlycosmos.clear();invalidate();if(glycosmosPageData)renderGlycosmos();else renderGlycosmosSelection();};
$('glycosmos-select-page').onclick=()=>{for(const g of glycosmosPageData.glycans)chooseGlycosmos(g,true);renderGlycosmos();};
$('glycosmos-clear').onclick=()=>{allGlycosmos=false;selectedGlycosmos.clear();excludedGlycosmos.clear();invalidate();if(glycosmosPageData)renderGlycosmos();else renderGlycosmosSelection();};
loadGlycosmos();

$('identify-glycan-adducts').onclick = async () => {
  const action=$('identify-glycan-adducts'), data=searchInputs(), signature=JSON.stringify(data);
  action.disabled=true;lastGlycanAdducts=null;$('glycan-adduct-results').hidden=true;
  status('Comparing glycan-only and adduct-associated signals…');
  try {
    const result=await api('/api/glycan-adducts',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)});
    if (JSON.stringify(searchInputs())!==signature) {status('Inputs changed. Identify glycan adducts again to update the calculation.');return;}
    lastGlycanAdducts=result;renderGlycanAdducts();status('Best glycan signal assignment complete.');
    $('glycan-adduct-results').scrollIntoView({behavior:'smooth',block:'start'});
  } catch(e) {status(e.message,true);}
  finally {action.disabled=false;}
};
function renderGlycanAdducts() {
  const data=lastGlycanAdducts, box=$('glycan-adduct-list');box.replaceChildren();
  box.append(el('p',`${data.unique_masses.toLocaleString()} selected free glycan masses compared. One best assignment is reported per m/z and charge, using the smallest absolute m/z error within tolerance (or the closest candidate if none fit). Glycan-only means normal protonation/deprotonation without an additional adduct.`, 'hint'));
  if (data.glycosmos_coverage) {
    const c=data.glycosmos_coverage;
    box.append(el('p',`GlyCosmos: ${c.selected.toLocaleString()} selected; ${c.supported.toLocaleString()} supported; ${c.skipped.toLocaleString()} skipped without a supported mass.`, 'hint'));
  }
  for (const result of data.results) {
    const block=el('article',undefined,'glycan-adduct-measurement'), candidate=result.best_assignment;
    const category=candidate.assignment_type==='glycan_only'?'Glycan-only':`Adduct-associated glycan · ${candidate.adduct}`;
    block.append(el('h3',`m/z ${result.mz} · charge ${result.charge}${data.ion_mode==='positive'?'+':'−'}`));
    const summary=el('p',`${candidate.within?'Best assignment':'Closest assignment (outside tolerance)'}: ${category}`,'note assignment-summary');
    summary.dataset.assignmentType=candidate.assignment_type;block.append(summary);
    if (!candidate.within) block.append(el('p','No tested explanation fits this measurement within your tolerance. This is the closest candidate, not an identification.','hint'));
    if (result.tied_assignment_count>1) block.append(el('p',`Ambiguous: ${result.tied_assignment_count} explanations tie within 0.000000001 m/z. One representative is shown; glycan-only is preferred on a tie. Mass alone cannot distinguish these explanations.`, 'note assignment-tie'));
    block.append(el('p',`Calculated MW = (${result.mz} × ${result.charge}) − ${result.charge} = ${fixed(result.calculated_mw)} Da`, 'equation'));
    block.append(el('p','This is the requested 1 Da-per-charge approximation.','hint'));
    block.append(el('h3',`Assigned glycan MW: ${fixed(candidate.glycan_mw)} Da`));
    block.append(el('p',`Mass difference = ${fixed(candidate.glycan_mw)} − ${fixed(result.calculated_mw)} = ${signed(candidate.mass_difference)} Da`, 'equation'));
    const references=el('div',undefined,'adduct-glycan-references');
    for (const id of candidate.glycan_info.ids) {
      if (/^G\d{5}[A-Z]{2}$/.test(id)) references.append(glycosmosReference(id,data));
      else {
        const glycan=catalog?.glycans.find(g=>g.id===id);
        if (glycan) {
          const ref=el('div',undefined,'glycosmos-result-reference');
          ref.append(el('p',glycan.name),glycanDiagram(glycan,true),button('Structure & subunit names ↗',()=>openStructure(glycan),'structure-button'));references.append(ref);
        }
      }
    }
    if (candidate.glycan_info.ids.length) {
      const details=el('details');details.open=true;
      details.append(el('summary',`Glycan references and structures (${candidate.glycan_info.ids.length} shown; up to 20 per mass)`),references);
      details.append(el('p','These references share the assigned mass; the structure is not uniquely identified by mass.','hint'));block.append(details);
    } else block.append(el('p','Custom glycan mass; no reference structure supplied.','hint'));
    block.append(el('p',`Ion form: ${candidate.ion} · multimer count n = ${candidate.multimer} · ${candidate.within?'Within':'Outside'} tolerance (${data.tolerance} ${data.unit}).`, 'hint'));
    block.append(el('p',`Expected difference = z − ion shift − (n − 1) × glycan MW = ${result.charge} − (${signed(candidate.mass_shift)}) − (${candidate.multimer} − 1) × ${fixed(candidate.glycan_mw)} = ${signed(candidate.expected_difference)} Da`, 'equation'));
    block.append(el('p',`Difference error = ${signed(candidate.mass_difference)} − (${signed(candidate.expected_difference)}) = ${signed(candidate.difference_error)} Da. Calculated m/z = ${fixed(candidate.predicted_mz)}.`, 'equation'));
    block.append(el('p',`Absolute m/z error: ${fixed(Math.abs(candidate.error_mz))} · signed ppm: ${signed(candidate.error_ppm)}.`, 'hint'));
    box.append(block);
  }
  $('glycan-adduct-results').hidden=false;
}
$('download-glycan-adducts').onclick=()=>{
  const data=lastGlycanAdducts;if(!data)return;
  const columns=['observed_mz','charge','assignment_type','tied_assignment_count','calculated_mw','assigned_glycan_mw','glycan_minus_calculated_mw','glycan_ids','glycan_names','ion','adduct','adduct_id','multimer','mass_shift','expected_difference','difference_error','predicted_mz','error_mz','error_ppm','error_da','within','ion_mode','tolerance','unit','formula','glycan_mass_basis','glycosmos_selection','selected_adducts'];
  const rows=[columns];
  for (const r of data.results) {
    const c=r.best_assignment;
    rows.push([r.mz,r.charge,c.assignment_type,r.tied_assignment_count,r.calculated_mw,c.glycan_mw,c.mass_difference,c.glycan_info.ids.join('; '),c.glycan_info.names.join('; '),...['ion','adduct','adduct_id','multimer','mass_shift','expected_difference','difference_error','predicted_mz','error_mz','error_ppm','error_da','within'].map(k=>c[k]??''),data.ion_mode,data.tolerance,data.unit,data.formula,data.glycan_mass_basis,JSON.stringify(data.glycosmos_selection),JSON.stringify(data.selected_adducts??'proton-only')]);
  }
  const quote=value=>'"'+String(value??'').replaceAll('"','""')+'"';
  const blob=new Blob([rows.map(r=>r.map(quote).join(',')).join('\r\n')],{type:'text/csv;charset=utf-8'});
  const url=URL.createObjectURL(blob),a=el('a');a.href=url;a.download='glypeptid-glycan-assignments.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
};
