// Optional browser smoke test: Node 22+ and Chrome running with
// --headless --remote-debugging-port=9223 --user-data-dir=/tmp/glycotag-browser-profile
// Start the Python app on port 8000 first. No npm packages required.
import assert from 'node:assert/strict';
import {writeFile} from 'node:fs/promises';
const pages = await (await fetch('http://127.0.0.1:9223/json/list')).json();
const target = pages.find(p=>p.type==='page');
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
  await send('Page.navigate',{url:process.env.GLYCOTAG_URL || 'http://localhost:8000/'});
  await waitFor('typeof catalog !== "undefined" && catalog !== null');
  assert.equal(await evaluate('document.querySelectorAll(".glycan-card").length'),21);
  await evaluate('$("accession").value="P02787";$("import").click()');
  await waitFor('!$("import").disabled && $("sites").value.includes("630,N")');
  assert.equal(await evaluate('$("sites").value'),'51,O\n432,N\n491,N\n630,N');
  assert.equal(await evaluate('document.querySelectorAll(".site-chip").length'),4);
  assert.match(await evaluate('$("import-info").textContent'),/atypical/);
  console.log('PASS: live UniProt import loads sequence, O/N sites and evidence.');
  await evaluate('document.querySelector(".glycan-group").open=true;document.querySelector(".glycan-card .structure-button").click()');
  assert.equal(await evaluate('$("structure-dialog").open'),true);
  assert.match(await evaluate('$("structure-content").textContent'),/N-acetylgalactosamine/);
  await screenshot('/tmp/glycotag-structure.png');
  await evaluate('$("close-structure").click();openStructure(catalog.glycans.find(g=>g.id==="n_a3g3s3"))');
  assert.match(await evaluate('$("structure-content").textContent'),/Mannose/);
  assert.equal(await evaluate('document.querySelectorAll("#structure-content .sugar-label").length'),14);
  await screenshot('/tmp/glycotag-complex-structure.png');
  await evaluate('$("close-structure").click();document.querySelector("input[value=o_core1]").click();$("mode").value="mono_attached";$("mode").dispatchEvent(new Event("change",{bubbles:true}))');
  assert.equal(await evaluate('document.querySelector("input[value=o_core1]").checked'),true);
  assert.match(await evaluate('document.querySelector("input[value=o_core1]").closest(".glycan-card").textContent'),/365\.13219594142/);
  console.log('PASS: reference structures, subunit names and mass convention conversion.');
  await evaluate(`
    function fill(id,value){$(id).value=value;$(id).dispatchEvent(new Event('input',{bubbles:true}));}
    fill('sequence','S');fill('sites','1,O');fill('glycans_O','');fill('glycans_N','');
    const m=catalog.glycans.find(g=>g.id==='o_core1').attached_mass+87.03202840472+18.01056468403;
    fill('observations',String(m/2-1.007276466621)+',2\\n1.1,1');
    $('form').requestSubmit();
  `);
  await waitFor('!$("results").hidden && !$("run").disabled');
  assert.equal(await evaluate('last.results.length'),2);
  assert.equal(await evaluate('last.results[0].match_count'),1);
  assert.equal(await evaluate('last.results[1].match_count'),0);
  assert.ok(Math.abs(await evaluate('last.results[0].rows[0].error_da'))<1e-8);
  assert.match(await evaluate('$("calculation-body").textContent'),/Serine/);
  assert.match(await evaluate('$("calculation-body").textContent'),/87\.03202840472/);
  assert.match(await evaluate('$("calculation-body").textContent'),/18\.01056468403/);
  assert.match(await evaluate('$("calculation-body").textContent'),/terminal H₂O/);
  assert.match(await evaluate('$("calculation-body").textContent'),/Calculated − target/);
  await evaluate('document.querySelectorAll("#result-list .result-block")[1].querySelector("button").click()');
  assert.equal(await evaluate('$("calc-measurement").value'),'1');
  assert.match(await evaluate('$("calculation-body").textContent'),/outside tolerance/);
  await evaluate('selectCalculation(0,0,true)');
  await new Promise(r=>setTimeout(r,500));
  await screenshot('/tmp/glycotag-calculation.png');
  console.log('PASS: batch search, exact match, no-match state and per-candidate calculation.');
  await evaluate('window.savedBlob=null;const originalURL=URL.createObjectURL;URL.createObjectURL=b=>{savedBlob=b;return originalURL(b)};$("download").click()');
  const csv=await evaluate('savedBlob.text()');
  assert.match(csv,/aa_water_loss/);assert.match(csv,/o_core1/);assert.match(csv,/mono_attached/);
  console.log('PASS: CSV contains calculation totals and glycan selection.');
  const batch = await evaluate("analyzeBatch({...inputs, observations:Array(21).fill(inputs.observations.split('\\n')[0]).join('\\n')}).then(r=>({count:r.results.length,matches:r.results.every(x=>x.match_count===1)}))");
  assert.equal(batch.count,21); assert.equal(batch.matches,true);
  console.log('PASS: large input is split into hosted batches and all 21 results are retained.');
  await evaluate(`const dt=new DataTransfer();dt.items.add(new File(['mass\\n700.123\\n800.456'], 'o-masses.csv',{type:'text/csv'}));$('file_O').files=dt.files;$('file_O').dispatchEvent(new Event('change',{bubbles:true}));`);
  await waitFor('$("glycans_O").value.includes("800.456")');
  assert.equal(await evaluate('$("results").hidden'),true);
  assert.equal(await evaluate('$("calculation").hidden'),true);
  assert.equal(await evaluate('$("glycans_O").value'),'700.123, 800.456');
  await evaluate(`const bad=new DataTransfer();bad.items.add(new File(['mass\\nNaN'], 'bad.csv',{type:'text/csv'}));$('file_O').files=bad.files;$('file_O').dispatchEvent(new Event('change',{bubbles:true}));`);
  await waitFor('$("mass-import-status").className==="error"');
  assert.equal(await evaluate('$("glycans_O").value'),'700.123, 800.456');
  console.log('PASS: custom CSV import, invalid-file rejection and stale result clearing.');
  await evaluate(`const realFetch=window.fetch;window.fetch=(url,...args)=>String(url).includes('/api/uniprot/')?Promise.resolve(new Response(JSON.stringify({sequence:'ASNT',accession:'TEST00',name:'Unannotated test protein',sites:[],sites_text:'',skipped:[],warnings:[],source:'https://www.uniprot.org/'}),{headers:{'Content-Type':'application/json'}})):realFetch(url,...args);$('accession').value='TEST00';$('import').click()`);
  await waitFor('!$("import").disabled');
  assert.equal(await evaluate('$("sites").value'),'');
  assert.match(await evaluate('$("import-info").textContent'),/does not mean/);
  console.log('PASS: empty annotation result clears previous sites and explains missing coverage.');
  await evaluate('$("example").click()');await waitFor('$("sites").value.includes("630,N")');
  await send('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});
  await evaluate('window.scrollTo(0,0)');await new Promise(r=>setTimeout(r,300));
  assert.ok(await evaluate('document.documentElement.scrollWidth <= window.innerWidth'), 'Mobile page has horizontal overflow');
  await screenshot('/tmp/glycotag-mobile.png');
  assert.deepEqual(errors,[]);
  console.log('PASS: mobile layout without page overflow; no browser JavaScript exceptions.');
} finally { ws.close(); }
