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
  await send('Runtime.enable');await send('Page.enable');
  await send('Emulation.setDeviceMetricsOverride',{width:1440,height:1050,deviceScaleFactor:1,mobile:false});
  await send('Page.navigate',{url:process.env.GLYPEPTID_URL || 'http://localhost:8000/'});
  await waitFor('typeof adductCatalog !== "undefined" && adductCatalog !== null && catalog !== null && glycosmosPageData !== null');
  await evaluate(`
    function fill(id,value){$(id).value=value;$(id).dispatchEvent(new Event('input',{bubbles:true}));}
    $('glycosmos-panel').open=true;
    $('glycosmos-query').value='G00039MO';loadGlycosmos();
  `);
  await waitFor('glycosmosPageData.glycans.some(g=>g.id==="G00039MO")');
  await evaluate(`
    document.querySelector('#glycosmos-library input[value=G00039MO]').click();
    document.querySelector('[data-mode=positive] .library-actions button').click();
    fill('tolerance','1');fill('observations','447.1585726,1');$('identify-glycan-adducts').click();
  `);
  await waitFor('lastGlycanAdducts !== null && !$("identify-glycan-adducts").disabled');
  assert.equal(await evaluate('$("sequence").value'),'');
  assert.equal(await evaluate('$("sites").value'),'');
  assert.equal(await evaluate('lastGlycanAdducts.results[0].calculated_mw'),446.1585726);
  assert.ok(Math.abs(await evaluate('lastGlycanAdducts.results[0].best_assignment.mass_difference')+21.989218)<1e-9);
  assert.equal(await evaluate('lastGlycanAdducts.results[0].best_assignment.adduct'),'[M+Na]+');
  assert.equal(await evaluate('lastGlycanAdducts.results[0].best_assignment.within'),true);
  assert.match(await evaluate('$("glycan-adduct-list").textContent'),/Calculated MW = \(447\.1585726 × 1\) − 1/);
  await evaluate('document.querySelector("#glycan-adduct-list .glycosmos-result-reference").scrollIntoView()');
  await waitFor('document.querySelector("#glycan-adduct-list img").naturalWidth>0');
  await evaluate('document.querySelector("#glycan-adduct-list .structure-button").click()');
  await waitFor('$("structure-dialog").open');
  assert.match(await evaluate('$("structure-content").textContent'),/galactosamine/);
  await evaluate('$("close-structure").click()');
  await evaluate('window.savedBlob=null;const originalURL=URL.createObjectURL;URL.createObjectURL=b=>{savedBlob=b;return originalURL(b)};$("download-glycan-adducts").click()');
  const csv=await evaluate('savedBlob.text()');
  assert.match(csv,/glycan_minus_calculated_mw/);assert.match(csv,/446\.1585726/);
  assert.match(csv,/G00039MO/);assert.match(csv,/\[M\+Na\]\+/);
  assert.match(csv,/expected_difference/);assert.match(csv,/assignment_type/);
  assert.equal(csv.split('\r\n').length,2);assert.match(csv,/adduct_associated/);
  await screenshot('/tmp/glypeptid-nearest-glycan-adducts.png');
  console.log('PASS: requested MW, signed nearest-glycan difference, sodium match, structures and CSV without protein input.');

  await evaluate(`fill('observations','447.1585726,1\\n235.0738953,2');$('identify-glycan-adducts').click()`);
  await waitFor('lastGlycanAdducts !== null && !$("identify-glycan-adducts").disabled');
  assert.equal(await evaluate('lastGlycanAdducts.results.length'),2);
  assert.equal(await evaluate('lastGlycanAdducts.results[1].best_assignment.adduct'),'[M+2Na]2+');
  assert.equal(await evaluate('lastGlycanAdducts.results[1].best_assignment.within'),true);
  await evaluate(`$('ion_mode').value='negative';$('ion_mode').dispatchEvent(new Event('change',{bubbles:true}));document.querySelector('[data-mode=negative] .library-actions button').click();fill('observations','1271.500787333379,1');$('identify-glycan-adducts').click()`);
  await waitFor('lastGlycanAdducts !== null && !$("identify-glycan-adducts").disabled');
  assert.equal(await evaluate('lastGlycanAdducts.results[0].best_assignment.adduct'),'[3M-H]-');
  assert.equal(await evaluate('lastGlycanAdducts.results[0].best_assignment.within'),true);
  await evaluate(`fill('observations','107,4');$('identify-glycan-adducts').click()`);
  await waitFor('lastGlycanAdducts !== null && !$("identify-glycan-adducts").disabled');
  assert.match(await evaluate('$("glycan-adduct-list").textContent'),/Closest assignment \(outside tolerance\): Glycan-only/);
  console.log('PASS: batch observations, multiply charged sodium, negative trimer and incompatible charge.');

  await evaluate(`$('glycosmos-clear').click();fill('glycans_O','1000,1022');$('ion_mode').value='positive';$('ion_mode').dispatchEvent(new Event('change',{bubbles:true}));fill('observations','1022.989218,1');$('identify-glycan-adducts').click()`);
  await waitFor('lastGlycanAdducts !== null && !$("identify-glycan-adducts").disabled');
  assert.equal(await evaluate('lastGlycanAdducts.results[0].best_assignment.glycan_mw'),1000);
  assert.match(await evaluate('$("glycan-adduct-list").textContent'),/Best assignment: Adduct-associated glycan/);
  assert.equal(await evaluate('document.querySelectorAll("#glycan-adduct-list .assignment-summary").length'),1);
  await evaluate(`fill('glycans_O','1000');fill('observations','1001.007276466621,1\\n1001.007276466621,1');$('identify-glycan-adducts').click()`);
  await waitFor('lastGlycanAdducts !== null && !$("identify-glycan-adducts").disabled');
  assert.equal(await evaluate('lastGlycanAdducts.results.length'),1);
  assert.equal(await evaluate('lastGlycanAdducts.results[0].best_assignment.assignment_type'),'glycan_only');
  assert.equal(await evaluate('lastGlycanAdducts.results[0].best_assignment.adduct'),null);
  assert.match(await evaluate('document.querySelector(".assignment-summary").textContent'),/Best assignment: Glycan-only/);
  await evaluate('$("download-glycan-adducts").click()');
  const glycanOnlyCsv=await evaluate('savedBlob.text()');
  assert.match(glycanOnlyCsv,/glycan_only/);assert.equal(glycanOnlyCsv.split('\r\n').length,2);
  await evaluate(`fill('glycans_O','1000,'+String(1022.989218-1.007276466621));fill('observations','1022.989218,1');$('identify-glycan-adducts').click()`);
  await waitFor('lastGlycanAdducts !== null && !$("identify-glycan-adducts").disabled');
  assert.equal(await evaluate('lastGlycanAdducts.results[0].tied_assignment_count'),2);
  assert.equal(await evaluate('document.querySelectorAll(".assignment-summary").length'),1);
  assert.match(await evaluate('document.querySelector(".assignment-tie").textContent'),/Ambiguous: 2 explanations/);
  assert.match(await evaluate('document.querySelector(".assignment-summary").textContent'),/Glycan-only/);
  await screenshot('/tmp/glypeptid-single-assignment-tie.png');
  await send('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});
  assert.ok(await evaluate('document.documentElement.scrollWidth <= window.innerWidth'),'Mobile overflow');
  await evaluate(`fill('observations','1,1')`);
  assert.equal(await evaluate('lastGlycanAdducts'),null);
  assert.equal(await evaluate('$("glycan-adduct-results").hidden'),true);
  await evaluate('$("identify-glycan-adducts").click()');
  await waitFor('!$("identify-glycan-adducts").disabled');
  assert.match(await evaluate('$("status").textContent'),/must give a positive finite MW/);
  assert.deepEqual(errors,[]);
  console.log('PASS: competing explanations, one CSV row, duplicate measurements, explicit ties, stale-result clearing, validation and mobile layout.');
} finally {ws.close();await fetch(debugURL+'/json/close/'+target.id);}
