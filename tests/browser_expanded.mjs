// Optional browser smoke test: Node 22+ and Chrome running with
// --headless --remote-debugging-port=9223 --user-data-dir=/tmp/glypeptid-browser-profile
// Start the Python app on port 8000 first. No npm packages required.
import assert from 'node:assert/strict';
import {writeFile} from 'node:fs/promises';
const debugURL = process.env.CHROME_DEBUG_URL || 'http://127.0.0.1:9223';
const target = await (await fetch(debugURL + '/json/new?about:blank', {method:'PUT'})).json();
const ws = new WebSocket(target.webSocketDebuggerUrl);
await new Promise((resolve,reject)=>{ws.onopen=resolve;ws.onerror=reject;});
let next=0; const pending=new Map(), errors=[];
ws.onmessage=e=>{const data=JSON.parse(e.data);if(data.id){const p=pending.get(data.id);pending.delete(data.id);if(data.error)p.reject(Error(data.error.message));else p.resolve(data.result);}else if(data.method==='Runtime.exceptionThrown')errors.push(data.params.exceptionDetails.text);};
function send(method,params={}){return new Promise((resolve,reject)=>{const id=++next;pending.set(id,{resolve,reject});ws.send(JSON.stringify({id,method,params}));});}
async function evaluate(expression){const r=await send('Runtime.evaluate',{expression,awaitPromise:true,returnByValue:true});if(r.exceptionDetails)throw Error(JSON.stringify(r.exceptionDetails));return r.result.value;}
async function waitFor(expression){const until=Date.now()+30000;while(Date.now()<until){if(await evaluate(expression))return;await new Promise(r=>setTimeout(r,100));}throw Error('Timed out: '+expression);}
async function screenshot(path){const shot=await send('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});await writeFile(path,Buffer.from(shot.data,'base64'));}
try {
  await send('Runtime.enable'); await send('Page.enable');
  await send('Emulation.setDeviceMetricsOverride',{width:1440,height:1050,deviceScaleFactor:1,mobile:false});
  await send('Page.navigate',{url:process.env.GLYPEPTID_URL || 'http://localhost:8000/'});
  await waitFor('typeof adductCatalog !== "undefined" && adductCatalog !== null && catalog !== null');
  assert.equal(await evaluate('$("adduct-library").hidden'),false);
  assert.equal(await evaluate('document.querySelectorAll(".adduct-group[open]").length'),2);
  assert.equal(await evaluate('document.querySelectorAll(".adduct-choice").length'),46);
  await evaluate(`document.querySelector('[data-mode="positive"] .library-actions button').click()`);
  assert.equal(await evaluate('document.querySelectorAll("[data-mode=positive] input:checked").length'),31);
  await evaluate(`document.querySelector('[data-mode="negative"] .library-actions button').click()`);
  assert.equal(await evaluate('document.querySelectorAll("[data-mode=negative] input:checked").length'),15);
  await evaluate(`document.querySelector('[data-mode="positive"] input').click()`);
  assert.equal(await evaluate('document.querySelectorAll("[data-mode=positive] input:checked").length'),30);
  await evaluate(`document.querySelector('[data-mode="positive"] .library-actions button').click()`);
  assert.equal(await evaluate('document.querySelectorAll("[data-mode=positive] input:checked").length'),31);
  console.log('PASS: all 46 adducts, both select-all sections and individual deselection.');

  assert.equal(await evaluate('$("glycosmos-panel").tagName'),'DETAILS');
  assert.equal(await evaluate('$("glycosmos-panel").open'),false);
  await evaluate('document.querySelector("#glycosmos-panel > summary").click()');
  assert.equal(await evaluate('$("glycosmos-panel").open'),true);
  await waitFor('glycosmosPageData !== null');
  assert.equal(await evaluate('glycosmosPageData.total'),257024);
  assert.equal(await evaluate('document.querySelectorAll(".glycosmos-card").length'),30);
  await evaluate('loadGlycosmos(2)');
  assert.equal(await evaluate('glycosmosPageData.page'),2);
  await evaluate(`$('glycosmos-query').value='G00039MO';$('glycosmos-supported').checked=true;loadGlycosmos()`);
  assert.ok(await evaluate('glycosmosPageData.glycans.some(g=>g.id==="G00039MO")'));
  await evaluate(`document.querySelector('#glycosmos-library input[value=G00039MO]').click()`);
  assert.equal(await evaluate('selectedGlycosmos.size'),1);
  assert.equal(await evaluate('selectedGlycosmos.has("G00039MO")'),true);
  await evaluate(`openGlycosmos(glycosmosPageData.glycans.find(g=>g.id==='G00039MO'))`);
  await waitFor('document.querySelector("#structure-content img").complete && document.querySelector("#structure-content img").naturalWidth>0');
  assert.match(await evaluate('$("structure-content").textContent'),/galactosamine/);
  await screenshot('/tmp/glypeptid-glycosmos-structure.png');
  await evaluate('$("close-structure").click()');
  await evaluate(`$('glycosmos-query').value='G98562EA';$('glycosmos-supported').checked=false;loadGlycosmos()`);
  assert.equal(await evaluate('glycosmosPageData.glycans[0].free_mass'),null);
  assert.match(await evaluate('$("glycosmos-library").textContent'),/No supported ExPASy mass/);
  assert.equal(await evaluate('document.querySelectorAll(".glycosmos-card input[type=checkbox]").length'),1);
  assert.equal(await evaluate('selectedGlycosmos.size'),1);
  console.log('PASS: full GlyCosmos library, pagination, search, selection persistence and unavailable masses.');

  await evaluate(`
    function fill(id,value){$(id).value=value;$(id).dispatchEvent(new Event('input',{bubbles:true}));}
    fill('sequence','S');fill('sites','1,O');fill('max_length','1');
    fill('observations',String(87.03202840472+424.1693546+22.989218)+',1');
    $('form').requestSubmit();
  `);
  await waitFor('last !== null && !$("run").disabled');
  assert.equal(await evaluate('last.results[0].rows[0].adduct'),'[M+Na]+');
  assert.ok(Math.abs(await evaluate('last.results[0].rows[0].absolute_delta_mz'))<1e-9);
  assert.match(await evaluate('$("calculation-body").textContent'),/G00039MO/);
  assert.match(await evaluate('$("calculation-body").textContent'),/22\.989218/);
  assert.match(await evaluate('$("calculation-body").textContent'),/511\.201383/);
  await evaluate('document.querySelector("#result-list .glycosmos-result-reference").scrollIntoView()');
  await waitFor('document.querySelector("#result-list .glycosmos-result-reference img").naturalWidth > 0');
  assert.equal(await evaluate('document.querySelector("#result-list .glycosmos-result-reference > a").textContent'),'G00039MO');
  assert.equal(await evaluate('document.querySelector("#calculation-body .glycosmos-result-reference img").getAttribute("src")'),'/api/glycosmos/image/G00039MO');
  // Clear the browse cache to exercise the exact-record endpoint from results.
  await evaluate('glycosmosRecords.clear();document.querySelector("#result-list .glycosmos-result-reference button").click()');
  await waitFor('$("structure-dialog").open');
  await waitFor('document.querySelector("#structure-content img").naturalWidth > 0');
  assert.match(await evaluate('$("structure-content").textContent'),/O-linked/);
  assert.match(await evaluate('$("structure-content").textContent'),/galactosamine/);
  await evaluate('$("close-structure").click()');
  await evaluate('window.savedBlob=null;const originalURL=URL.createObjectURL;URL.createObjectURL=b=>{savedBlob=b;return originalURL(b)};$("download").click()');
  const csv=await evaluate('savedBlob.text()');
  assert.match(csv,/adduct_mass_shift/); assert.match(csv,/multimer/); assert.match(csv,/G00039MO/);
  assert.match(csv,/\[M\+Na\]\+/); assert.match(csv,/glycosmos_selection/);
  await screenshot('/tmp/glypeptid-expanded-calculation.png');
  await evaluate(`$('glycan_mass_type').value='attached';$('glycan_mass_type').dispatchEvent(new Event('change',{bubbles:true}));$('form').requestSubmit()`);
  await waitFor('last !== null && !$("run").disabled');
  assert.ok(Math.abs(await evaluate('last.results[0].rows[0].absolute_delta_mz'))<1e-9);
  assert.match(await evaluate('$("calculation-body").textContent'),/406\.1587899/);
  console.log('PASS: GlyCosmos sodium match, free/attached equivalence, adduct-specific calculation and CSV.');

  await evaluate(`$('ion_mode').value='negative';$('ion_mode').dispatchEvent(new Event('change',{bubbles:true}));fill('observations',String(3*(87.03202840472+424.1693546)-1.007276466621)+',1');$('form').requestSubmit()`);
  await waitFor('last !== null && !$("run").disabled');
  assert.equal(await evaluate('last.results[0].rows[0].adduct'),'[3M-H]-');
  assert.equal(await evaluate('last.results[0].rows[0].multimer'),3);
  assert.ok(Math.abs(await evaluate('last.results[0].rows[0].absolute_delta_mz'))<1e-9);
  assert.equal(await evaluate('document.querySelectorAll("[data-mode=positive] input:checked").length'),31);
  assert.match(await evaluate('$("calculation-body").textContent'),/multimer count n = 3/);
  await evaluate(`fill('observations','400,4');$('form').requestSubmit()`);
  await waitFor('last !== null && !$("run").disabled');
  assert.equal(await evaluate('last.results[0].rows.length'),0);
  assert.match(await evaluate('$("calculation-body").textContent'),/No selected adduct/);
  await evaluate(`document.querySelector('[data-mode="negative"] .library-actions button:nth-child(2)').click();$('form').requestSubmit()`);
  await waitFor('!$("run").disabled');
  assert.equal(await evaluate('last'),null);
  assert.match(await evaluate('$("status").textContent'),/Select at least one negative/);
  console.log('PASS: negative trimer, charge mismatch and empty-selection validation.');

  await evaluate(`$('glycosmos-all').click()`);
  assert.equal(await evaluate('$("glycosmos-count").textContent'),'257,024 SELECTED');
  // Select-all spans the entire catalogue even while a one-record filter is active.
  await evaluate(`document.querySelector('#glycosmos-library input').click()`);
  assert.equal(await evaluate('$("glycosmos-count").textContent'),'257,023 SELECTED');
  await evaluate('loadGlycosmos(1)');
  assert.equal(await evaluate('document.querySelector("#glycosmos-library input").checked'),false);
  await evaluate(`$('glycosmos-all').click();$('adduct-mode').value='proton';$('adduct-mode').dispatchEvent(new Event('change',{bubbles:true}));$('form').requestSubmit()`);
  await waitFor('last !== null && !$("run").disabled');
  assert.equal(await evaluate('last.glycosmos_coverage.selected'),257024);
  assert.equal(await evaluate('last.glycosmos_coverage.supported'),116089);
  assert.match(await evaluate('$("result-list").textContent'),/140,935 skipped/);
  await evaluate('$("download").click()');
  const allCsv=await evaluate('savedBlob.text()');
  assert.match(allCsv,/glycosmos_coverage/);assert.match(allCsv,/257024/);
  await evaluate('$("glycosmos-clear").click()');
  assert.equal(await evaluate('$("glycosmos-count").textContent'),'0 SELECTED');
  assert.equal(await evaluate('last'),null);
  console.log('PASS: all-record selection, exclusions across pages, full-library search and coverage export.');

  await evaluate(`$('glycosmos-query').value='';$('glycosmos-supported').checked=false;fill('sequence','NS');fill('sites','1,N\\n2,O');`);
  for (const [kind,total] of [['N',22667],['O',8850],['both',31517]]) {
    await evaluate(`$('glycosmos-linkage').value='${kind}';$('glycosmos-linkage').dispatchEvent(new Event('change',{bubbles:true}));`);
    await waitFor(`glycosmosPageData?.linkage === '${kind}' && !$('glycosmos-all').disabled`);
    assert.equal(await evaluate('glycosmosPageData.total'),total);
    assert.equal(await evaluate('selectedGlycosmos.size'),0);
    assert.equal(await evaluate('allGlycosmos'),false);
    assert.ok(await evaluate(`glycosmosPageData.glycans.every(g=>g.linkage_types.some(k=>'${kind}'==='both'||k==='${kind}'))`));
    await evaluate(`$('glycosmos-all').click()`);
    assert.equal(await evaluate('$("glycosmos-count").textContent'),`${total.toLocaleString('en-US')} SELECTED`);
    await evaluate(`$('form').requestSubmit()`);
    await waitFor('last !== null && !$("run").disabled');
    assert.equal(await evaluate('last.glycosmos_coverage.selected'),total);
    assert.equal(await evaluate('last.glycosmos_coverage.linkage'),kind);
    await evaluate(`document.querySelector('#glycosmos-library input').click()`);
    assert.equal(await evaluate('$("glycosmos-count").textContent'),`${(total-1).toLocaleString('en-US')} SELECTED`);
    await evaluate('loadGlycosmos(2)');
    assert.equal(await evaluate('glycosmosPageData.page'),2);
    await evaluate('loadGlycosmos(1)');
    assert.equal(await evaluate('document.querySelector("#glycosmos-library input").checked'),false);
  }
  const batchImages = await evaluate(`(async()=>{
    const data=searchInputs();
    const high=await api('/api/analyze',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({...data,observations:'4000,1'})});
    const combined=await analyzeBatch({...data,observations:[...Array(10).fill('400,1'),'4000,1'].join('\\n')});
    const firstIds=new Set(combined.results[0].rows.flatMap(r=>r.glycan_info.ids));
    return {count:combined.results.length, distinct:high.glycosmos_image_ids.some(id=>!firstIds.has(id)), preserved:high.glycosmos_image_ids.every(id=>combined.glycosmos_image_ids.includes(id))};
  })()`);
  assert.equal(batchImages.count,11);
  assert.equal(batchImages.distinct,true);
  assert.equal(batchImages.preserved,true);
  await evaluate(`$('glycosmos-linkage').value='N';$('glycosmos-linkage').dispatchEvent(new Event('change',{bubbles:true}));$('glycosmos-query').value='G00039MO';loadGlycosmos()`);
  assert.equal(await evaluate('glycosmosPageData.total'),0);
  await evaluate(`$('glycosmos-linkage').value='all';$('glycosmos-linkage').dispatchEvent(new Event('change',{bubbles:true}));$('glycosmos-query').value='';loadGlycosmos()`);
  assert.equal(await evaluate('glycosmosPageData.total'),257024);
  console.log('PASS: N/O/both source filters, scoped select-all, exclusions, pagination and selection reset.');

  await evaluate(`$('glycosmos-query').value='';loadGlycosmos()`);
  await evaluate(`$('glycosmos-panel').scrollIntoView();`);
  await screenshot('/tmp/glypeptid-expanded-library.png');
  await send('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});
  assert.ok(await evaluate('document.documentElement.scrollWidth <= window.innerWidth'),'Mobile overflow');
  await evaluate(`$('adduct-panel').scrollIntoView()`);
  await screenshot('/tmp/glypeptid-adducts-mobile.png');
  assert.deepEqual(errors,[]);
  console.log('PASS: expanded mobile layout; no JavaScript errors.');
} finally { ws.close(); await fetch(debugURL + '/json/close/' + target.id); }
