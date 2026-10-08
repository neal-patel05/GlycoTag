// Local-only browser checks. Start app.py --port 8001 and headless Chrome on 9223.
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

async function upload(name, text) {
  await evaluate(`(() => {const dt=new DataTransfer();dt.items.add(new File([${JSON.stringify(text)}],${JSON.stringify(name)},{type:'text/plain'}));$('peak-file').files=dt.files;$('peak-file').dispatchEvent(new Event('change'));})()`);
  await waitFor('!$("peak-status").textContent.startsWith("Reading")');
}
try {
  await send('Runtime.enable'); await send('Page.enable');
  await send('Emulation.setDeviceMetricsOverride',{width:1440,height:1100,deviceScaleFactor:1,mobile:false});
  await send('Page.navigate',{url:'http://localhost:8001/'});
  await waitFor('typeof catalog !== "undefined" && catalog !== null');
  assert.match(await evaluate('document.title'), /GlyPeptID/);
  assert.equal(await evaluate('document.querySelector(".brand img").naturalWidth'),428);
  await evaluate(`const dt=new DataTransfer();dt.items.add(new File(['not a real RAW'], 'sample.raw'));$('raw-file').files=dt.files;$('raw-file').dispatchEvent(new Event('change'));`);
  assert.match(await evaluate('$("raw-status").textContent'),/not been uploaded or analyzed/);
  await upload('qual.tsv', 'Mass\tIntensity\tCharge\n1000.123456789\t900\t3\n700.25\t500\t2\n612.34\t100\t0\n');
  assert.equal(await evaluate('document.querySelectorAll("#peak-rows tr").length'),3);
  assert.match(await evaluate('$("peak-status").textContent'),/1 have unknown charge/);
  assert.equal(await evaluate('$("peak-add").disabled'),true);
  await evaluate(`$('peak-ion-mode').value='positive';$('peak-ion-mode').dispatchEvent(new Event('change'));$('peak-select').click();$('peak-add').click();`);
  assert.equal(await evaluate('$("observations").value'), '1000.123456789,3\n700.25,2');
  assert.equal(await evaluate('document.querySelectorAll("#peak-rows input:checked").length'),2);
  await evaluate(`$('peak-add').click()`);
  assert.match(await evaluate('$("peak-status").textContent'),/already in your measurements/);
  await screenshot('/tmp/glypeptid-desktop.png');
  await upload('negative.csv', 'm/z,Charge\n900.75,3-');
  assert.equal(await evaluate('$("peak-ion-mode").value'),'negative');
  await evaluate(`$('peak-select').click();$('peak-add').click()`);
  assert.match(await evaluate('$("peak-status").textContent'),/differs from your existing measurements/);
  assert.equal(await evaluate('$("observations").value'), '1000.123456789,3\n700.25,2');
  await evaluate(`$('observations').value='';$('peak-add').click()`);
  assert.equal(await evaluate('$("observations").value'),'900.75,3');
  assert.equal(await evaluate('$("ion_mode").value'),'negative');
  await upload('mixed.csv','m/z,Z\n600,2+\n601,3-\n602,4');
  assert.equal(await evaluate('$("peak-ion-mode").value'),'');
  await evaluate(`$('peak-ion-mode').value='positive';$('peak-ion-mode').dispatchEvent(new Event('change'));$('peak-select').click();`);
  assert.equal(await evaluate('document.querySelectorAll("#peak-rows input:checked").length'),2);
  await upload('many.csv','m/z,Charge\n'+Array.from({length:101},(_,i)=>`${500+i},2`).join('\n'));
  await evaluate(`$('peak-ion-mode').value='positive';$('peak-ion-mode').dispatchEvent(new Event('change'));$('peak-select').click();$('peak-next').click();$('peak-select').click();$('peak-next').click();$('peak-select').click()`);
  assert.match(await evaluate('$("peak-add").textContent'),/100 selected/);
  await evaluate(`$('ion_mode').value='positive';$('peak-add').click()`);
  assert.match(await evaluate('$("peak-status").textContent'),/exceed 100/);
  assert.equal(await evaluate('$("observations").value'),'900.75,3');
  await upload('invalid.csv', 'Mass,Intensity\n500,2');
  assert.match(await evaluate('$("peak-status").textContent'),/does not supply charge/);
  assert.equal(await evaluate('$("peak-review").hidden'),true);
  assert.equal(await evaluate('$("observations").value'),'900.75,3');
  const exactMz = (87.03202840472 + 674.2382) / 2 + 1.007276466621;
  await upload('match.csv',`m/z,Charge\n${exactMz},2+`);
  await evaluate(`$('observations').value='';$('peak-select').click();$('peak-add').click();$('sequence').value='S';$('sites').value='1,O';$('glycans_O').value='674.2382';$('form').requestSubmit()`);
  await waitFor('!$("results").hidden && !$("run").disabled');
  assert.equal(await evaluate('last.results[0].match_count'),1);
  assert.ok(Math.abs(await evaluate('last.results[0].rows[0].signed_delta_mz')) < 1e-10);
  await upload('qual.csv','m/z,Charge\n600.25,2');
  await send('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});
  assert.equal(await evaluate('document.documentElement.scrollWidth <= innerWidth'),true);
  await evaluate('window.scrollTo(0,0)');
  await screenshot('/tmp/glypeptid-mobile.png');
  assert.deepEqual(errors,[]);
  console.log('PASS: logo, honest RAW state, peak preview, missing charge, polarity handling, duplicates, pagination, limits, failed imports, and mobile layout.');
} finally { ws.close(); }
