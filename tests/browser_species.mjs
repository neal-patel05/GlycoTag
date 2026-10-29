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
  await waitFor('typeof glycosmosPageData !== "undefined" && glycosmosPageData !== null && catalog !== null && adductCatalog !== null');
  await evaluate(`
    function fill(id,value){$(id).value=value;$(id).dispatchEvent(new Event('input',{bubbles:true}));}
    function linkage(value){$('glycosmos-linkage').value=value;$('glycosmos-linkage').dispatchEvent(new Event('change',{bubbles:true}));}
    $('glycosmos-panel').open=true;
    fill('glycosmos-species','humans');
  `);
  await waitFor('glycosmosPageData?.species === "humans"');
  assert.equal(await evaluate('glycosmosPageData.total'),8016);
  assert.match(await evaluate('$("glycosmos-species-status").textContent'),/Homo sapiens \(9606\)/);
  assert.equal(await evaluate('$("glycosmos-species-options").firstChild.value'),'Homo sapiens');
  assert.equal(await evaluate('glycosmosPageData.glycans.every(g=>g.organisms.some(t=>t.id==="9606"))'),true);
  await evaluate('$("glycosmos-all").click()');
  assert.equal(await evaluate('$("glycosmos-count").textContent'),'8,016 SELECTED');
  await evaluate('document.querySelector(".glycosmos-card input").click()');
  assert.equal(await evaluate('$("glycosmos-count").textContent'),'8,015 SELECTED');
  await evaluate('loadGlycosmos(2)');
  assert.equal(await evaluate('glycosmosPageData.page'),2);
  assert.equal(await evaluate('excludedGlycosmos.size'),1);
  for (const [kind,count] of [['N',2425],['O',892],['both',3317]]) {
    await evaluate(`linkage('${kind}')`);
    await waitFor(`glycosmosPageData?.linkage === '${kind}'`);
    assert.equal(await evaluate('glycosmosPageData.total'),count);
    assert.equal(await evaluate('$("glycosmos-count").textContent'),'0 SELECTED');
    await evaluate('$("glycosmos-all").click()');
    assert.equal(await evaluate('$("glycosmos-count").textContent'),count.toLocaleString('en-US')+' SELECTED');
  }
  // Immediately invalidate selection and old cards, even before debounce finishes.
  const reset=await evaluate(`fill('glycosmos-species','bovine');[allGlycosmos,selectedGlycosmos.size,excludedGlycosmos.size,$('glycosmos-all').disabled,$('glycosmos-select-page').disabled,document.querySelectorAll('.glycosmos-card').length]`);
  assert.deepEqual(reset,[false,0,0,true,true,0]);
  await waitFor('glycosmosPageData?.species === "bovine"');
  assert.equal(await evaluate('glycosmosPageData.total'),950);
  await evaluate(`linkage('all')`);
  await waitFor('glycosmosPageData?.linkage === "all"');
  assert.equal(await evaluate('glycosmosPageData.total'),1755);
  assert.match(await evaluate('$("glycosmos-species-status").textContent'),/Bos taurus \(9913\)/);
  console.log('PASS: human/bovine aliases, suggestions, all/N/O/both filters, select-all, exclusions and pagination.');

  await evaluate(`$('glycosmos-query').value='G00039MO';$('glycosmos-supported').checked=true;loadGlycosmos()`);
  await evaluate('$("glycosmos-all").click()');
  assert.equal(await evaluate('$("glycosmos-count").textContent'),'1,755 SELECTED');
  await evaluate('$("glycosmos-clear").click();document.querySelector("input[value=G00039MO]").click()');
  assert.match(await evaluate('document.querySelector("input[value=G00039MO]").closest("article").textContent'),/Bos taurus/);
  await evaluate(`openGlycosmos(glycosmosRecords.get('G00039MO'))`);
  await waitFor('document.querySelector("#structure-content img").naturalWidth > 0');
  assert.match(await evaluate('$("structure-content").textContent'),/Escherichia coli O45:H2/);
  await evaluate('$("close-structure").click()');
  await evaluate(`fill('observations','425.176631066621,1');$('identify-glycan-adducts').click()`);
  await waitFor('lastGlycanAdducts !== null && !$("identify-glycan-adducts").disabled');
  assert.equal(await evaluate('lastGlycanAdducts.glycosmos_selection.species'),'bovine');
  assert.equal(await evaluate('lastGlycanAdducts.glycosmos_coverage.selected'),1);
  assert.deepEqual(await evaluate('lastGlycanAdducts.results[0].best_assignment.glycan_info.ids'),['G00039MO']);
  await evaluate(`window.savedBlob=null;const originalURL=URL.createObjectURL;URL.createObjectURL=b=>{savedBlob=b;return originalURL(b)};$('download-glycan-adducts').click()`);
  assert.match(await evaluate('savedBlob.text()'),/bovine/);
  await evaluate(`fill('sequence','S');fill('sites','1,O');fill('max_length','1');fill('observations','512.208659471341,1');$('form').requestSubmit()`);
  await waitFor('last !== null && !$("run").disabled');
  assert.equal(await evaluate('last.glycosmos_coverage.species_matches[0].name'),'Bos taurus');
  await evaluate('$("download").click()');
  assert.match(await evaluate('savedBlob.text()'),/bovine/);
  console.log('PASS: source organism labels, structures, both calculation modes and CSV species provenance.');

  // A delayed old response must never replace a more recent species scope.
  await evaluate(`
    $('glycosmos-query').value='';$('glycosmos-supported').checked=false;
    window.originalSpeciesFetch=window.fetch;window.speciesDelayStarted=false;
    window.fetch=async(...args)=>{
      const response=await originalSpeciesFetch(...args);
      if(String(args[0]).includes('species=human&') || String(args[0]).endsWith('species=human')) {
        speciesDelayStarted=true;await new Promise(r=>setTimeout(r,1200));
      }
      return response;
    };
    fill('glycosmos-species','human');loadGlycosmos();
  `);
  await waitFor('speciesDelayStarted');
  await evaluate(`fill('glycosmos-species','bovine');loadGlycosmos()`);
  await waitFor('glycosmosPageData?.species === "bovine"');
  await new Promise(r=>setTimeout(r,1400));
  assert.equal(await evaluate('glycosmosPageData.species'),'bovine');
  await evaluate('window.fetch=originalSpeciesFetch');
  await evaluate(`fill('glycosmos-species','not-an-organism')`);
  await waitFor('glycosmosPageData?.species === "not-an-organism"');
  assert.equal(await evaluate('glycosmosPageData.total'),0);
  assert.equal(await evaluate('$("glycosmos-all").disabled'),true);
  assert.match(await evaluate('$("glycosmos-species-status").textContent'),/No organism annotations/);
  assert.equal(await evaluate('last'),null);
  assert.equal(await evaluate('lastGlycanAdducts'),null);
  await evaluate(`fill('glycosmos-species','')`);
  await waitFor('glycosmosPageData?.species === ""');
  assert.equal(await evaluate('glycosmosPageData.total'),257024);
  await evaluate(`fill('glycosmos-species','bovine')`);
  await waitFor('glycosmosPageData?.species === "bovine"');
  await evaluate('$("glycosmos-species").scrollIntoView()');
  await screenshot('/tmp/glypeptid-species-desktop.png');
  await send('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});
  assert.ok(await evaluate('document.documentElement.scrollWidth <= window.innerWidth'),'Mobile overflow');
  await evaluate('$("glycosmos-species").scrollIntoView()');
  await screenshot('/tmp/glypeptid-species-mobile.png');
  assert.deepEqual(errors,[]);
  console.log('PASS: stale-response protection, no-match state, clearing species, mobile layout and no JavaScript errors.');
} finally {ws.close();await fetch(debugURL + '/json/close/' + target.id);}
