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
  await send('Page.navigate',{url:process.env.GLYPEPTID_URL || 'http://localhost:8000/'});
  await waitFor('typeof glycosmosPageData !== "undefined" && glycosmosPageData !== null && catalog !== null && adductCatalog !== null');
  await evaluate(`
    function fill(id,value){$(id).value=value;$(id).dispatchEvent(new Event('input',{bubbles:true}));}
    async function choose(id){$('glycosmos-query').value=id;await loadGlycosmos();document.querySelector('#glycosmos-library input[value='+id+']').click();}
    $('glycosmos-panel').open=true;
    fill('sequence','S');fill('sites','1,O');fill('max_length','1');
    fill('observations','998.367059471341,1');
    choose('G00026MO');
  `);
  await waitFor('selectedGlycosmos.has("G00026MO")');
  await evaluate('$("form").requestSubmit()');
  await waitFor('!$("run").disabled');
  assert.equal(await evaluate('last'),null);
  assert.equal(await evaluate('$("results").hidden'),true);
  assert.match(await evaluate('$("status").textContent'),/No linkage-compatible/);
  console.log('PASS: N-linked glycan at an O site is rejected, even at the exact expected mass.');

  await evaluate(`fill('sequence','NS');fill('sites','1,N\\n2,O');$('form').requestSubmit()`);
  await waitFor('last !== null && !$("run").disabled');
  assert.equal(await evaluate('last.results[0].rows.every(r=>r.type==="N")'),true);
  assert.match(await evaluate('$("result-list").textContent'),/Sites skipped.*O · 2/);
  await evaluate(`choose('G00039MO')`);
  await evaluate(`fill('sequence','S');fill('sites','1,O');$('form').requestSubmit()`);
  await waitFor('last !== null && !$("run").disabled');
  assert.deepEqual(await evaluate('last.results[0].rows.flatMap(r=>r.glycan_info.ids)'),['G00039MO']);
  assert.equal(await evaluate('last.glycosmos_coverage.incompatible'),1);
  assert.equal(await evaluate('document.querySelectorAll("#result-list a[href$=G00026MO]").length'),0);
  assert.match(await evaluate('$("result-list").textContent'),/1 additional records skipped/);
  await evaluate(`$('glycosmos-apply-to').value='N';$('glycosmos-apply-to').dispatchEvent(new Event('change',{bubbles:true}));$('form').requestSubmit()`);
  await waitFor('!$("run").disabled');
  assert.equal(await evaluate('last'),null);
  assert.match(await evaluate('$("status").textContent'),/No linkage-compatible/);
  console.log('PASS: mixed-site searches skip unavailable linkage; site setting cannot force N glycans onto O sites.');

  await evaluate(`$('glycosmos-clear').click();$('glycosmos-apply-to').value='both';choose('G00017IP')`);
  await evaluate(`choose('G04143WH')`);
  await evaluate(`fill('sequence','NS');fill('sites','1,N\\n2,O');fill('observations','1484.482659471341,1');$('form').requestSubmit()`);
  await waitFor('last !== null && !$("run").disabled');
  assert.equal(await evaluate('last.results[0].rows.every(r=>r.glycan_info.ids.length===1 && r.glycan_info.ids[0]===(r.type==="N"?"G00017IP":"G04143WH"))'),true);
  assert.equal(await evaluate('last.results[0].rows[0].type'),'O');
  assert.equal(await evaluate('last.results[0].rows[0].within'),true);
  await evaluate(`window.savedBlob=null;const originalURL=URL.createObjectURL;URL.createObjectURL=b=>{savedBlob=b;return originalURL(b)};$('download').click()`);
  const csv=await evaluate('savedBlob.text()');
  assert.match(csv,/linkage_enforced/);
  assert.match(csv,/G04143WH/);
  const csvRows=csv.split('\r\n').slice(1).filter(Boolean);
  for (const row of csvRows) {
    // The glycan_ids field precedes selection JSON (which legitimately lists both).
    const cells=row.match(/"(?:[^"]|"")*"/g).map(s=>s.slice(1,-1).replaceAll('""','"'));
    const header=csv.split('\r\n')[0].split(',').map(s=>s.slice(1,-1));
    const kind=cells[header.indexOf('type')], ids=cells[header.indexOf('glycan_ids')];
    assert.equal(ids,kind==='N'?'G00017IP':'G04143WH');
  }
  await evaluate('document.querySelector("#result-list .glycosmos-result-reference button").click()');
  await waitFor('$("structure-dialog").open');
  assert.match(await evaluate('$("structure-content").textContent'),/O-linked/);
  assert.match(await evaluate('$("structure-title").textContent'),/G04143WH/);
  await evaluate('$("close-structure").click()');
  assert.deepEqual(errors,[]);
  console.log('PASS: identical-mass N/O references stay separate in results, structure details and CSV.');
} finally {ws.close();await fetch(debugURL + '/json/close/' + target.id);}
