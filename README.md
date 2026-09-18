# GlycoTag

**Live website: https://glycotag.vercel.app** · [Vercel dashboard](https://vercel.com/neals-team/glycotag)

A Python website for ranking single-glycan, contiguous peptide candidates against observed m/z and charge. Python 3.9+; no packages required. The browser interface uses HTML/CSS/JavaScript, with all mass calculations in Python.

## Command reference

Run commands from the project folder:

```sh
cd "/Users/nealpatel/neal paper"
```

### Start locally

```sh
python3 app.py
```

Open **http://localhost:8000**. Leave that terminal running while using the local website. No Python package installation is needed.

### Stop locally

Press **Ctrl+C** in the terminal running the server (hold Control and press C; on a Mac, use Control, not Command). A `KeyboardInterrupt` message is normal. Closing the browser tab does not stop the server.

Stopping the local server does **not** stop the public website at https://glycotag.vercel.app.

### Restart after changes

Press **Ctrl+C**, then run:

```sh
python3 app.py
```

Python changes require restarting the server. For HTML, CSS, or JavaScript changes, refresh the browser; use **Command+Shift+R** on a Mac if old assets are cached. The local server does not automatically reload Python code.

### Use another port

```sh
python3 app.py --port 8001
```

Then open **http://localhost:8001**. Stop this server with **Ctrl+C** too.

### If port 8000 is already in use

On macOS, identify the process listening on the port:

```sh
lsof -nP -iTCP:8000 -sTCP:LISTEN
```

If it is the GlycoTag Python server you want to stop, note its `PID` and run `kill` with that number. For example, **only if the displayed PID is 12345**:

```sh
kill 12345
```

Then start the server again. Check the process before stopping it; alternatively, use port 8001 to leave the existing process running.

### Test the code

```sh
python3 -m unittest discover -s tests -v
node --check static/main.js
```

The second command checks JavaScript syntax and requires Node.js. See [Browser checks](#browser-checks) below for interactive functionality tests.

### Check that a server responds

With the local server running, use a second terminal:

```sh
curl --fail --silent --show-error http://localhost:8000/api/example
```

To check the public website:

```sh
curl --head https://glycotag.vercel.app
```

### Deploy updates to the public website

Requires Node.js/npm and access to the Vercel account. If you are not signed in:

```sh
npx vercel login
```

Publish the current files to the existing project:

```sh
npx vercel --prod --scope neals-team
```

Wait for a successful deployment. The public link remains **https://glycotag.vercel.app**. Local edits do not reach the public website until you deploy.

For a preview deployment that does not replace production:

```sh
npx vercel --scope neals-team
```

Open the preview URL printed by the CLI. Preview URLs may require Vercel authentication.

This folder is already linked to `neals-team/glycotag` through its local `.vercel/` metadata. Keep that folder locally. If working from a fresh copy and prompted, link to the **existing** `glycotag` project under **Neals Team**. More details: [DEPLOYMENT.md](DEPLOYMENT.md).

### GitHub and automatic deployment

This project was deployed directly with the Vercel CLI. **GitHub auto-deployment has not been set up.** Inspect the local Git state and configured remotes with:

```sh
git status --short
git diff
git remote -v
```

An empty `git remote -v` output means no remote is configured in this local repository. `git diff` shows tracked-file changes; use `git status` to also see new files. Connecting a GitHub repository to the Vercel project is a separate setup step. Until then, use the manual production deployment command above.

## Using the website

1. Click **Load serotransferrin**, paste one raw/FASTA protein sequence, or import a UniProt accession / entry ID (e.g. `P02787` or `TRFE_HUMAN`).
2. UniProt import fills in the sequence and its annotated N/O glycosylation sites. Review the descriptions and evidence, and edit or add sites as needed, one `position,O` or `position,N` per line. Positions are 1-based in the supplied sequence, including any signal peptide.
3. Enter measurements as `m/z,charge`, one per line, with no header. Charge is a positive integer.
4. Set tolerance, peptide length range, and calculation model. Browse the common glycan library and select entries to search, use **Select all glycans** to include all 21 O/N-linked references, or supply custom masses. **Clear selection** removes all library selections.
5. Run the search. Use **View calculation** for any candidate to see the amino acids, water corrections, glycan contribution, target mass, and signed difference. Download the displayed candidates and calculation totals as CSV.

The supplied 698-residue serotransferrin example contains S51 and N432/N491/N630. The complex/atypical annotations do not establish which glycan mass belongs to a site; selected N-linked masses are evaluated at all N sites.

## Automatic site import

The importer reads UniProt's JSON glycosylation annotations, preserving descriptions and evidence codes. It includes exact N-linked sites on N and O-linked sites on S/T, including atypical N sites. Evidence may be experimental or inferred; these are database annotations, not new predictions. Uncertain locations, ranges, other linkage types, and unsupported residues are listed as skipped. Missing annotations do not prove that a protein lacks glycosylation; enter known sites manually when needed.

Every new import replaces the prior sequence and sites. If an accession specifies an isoform whose sequence differs from the canonical annotated sequence, its FASTA sequence is imported without transferring canonical positions. The interface explains that isoform-specific sites must be entered manually. Editing a sequence or its sites marks the imported context as edited; changing search inputs clears outdated results.

## Glycan library and custom masses

The library contains 21 representative mammalian glycans: Tn, core 1/T antigen, sialyl-Tn, sialylated and disialylated core 1, core 2, disialylated extended core 2, core 3, O-GlcNAc, the trimannosyl N core, Man5, Man9, G0/G0F/G1F/G2/G2F, A2G2S1, A2G2S2, FA2G2S2, and A3G3S3. Entries have vector structure diagrams, linkage labels, named sugar subunits and counts, and reference links. Select **Structure & subunit names** to enlarge a diagram. These structures describe the library references; mass cannot uniquely identify a glycan or its isomer.

Library free glycan masses match [ExPASy GlycanMass](https://web.expasy.org/glycanmass/) with underivatized, monoisotopic settings. All 21 compositions were verified against the calculator on 2026-09-18; results are saved in `tests/fixtures/expasy-glycanmass.json`. It uses residue masses of Hex 162.0528, HexNAc 203.0794, Fuc 146.0579, and NeuAc 291.0954 Da, plus a free reducing end of 18.0105546 Da. Attached increments subtract the engine’s configured water (18.010564684 Da) from that free mass to preserve identical glycopeptide totals across models. Compact previews show sugars only; **Structure & subunit names** also shows the amino acid attachment. Library selections automatically supply a free glycan mass for `user` / `mono_free`, or an attached increment for `mono_attached`. Isobaric library entries share a mass search while retaining all their selected names and diagrams.

Custom fields start empty; enter, paste, or import masses as needed. Custom masses are used exactly as entered and **are not converted** when changing models. They are not assigned structures based on approximate agreement with the library. Library selections are additional to any custom values.

Paste comma-/whitespace-separated custom masses or use the O/N file inputs. TXT/CSV files contain only numeric mass values separated by commas, semicolons, or whitespace, with optional first header `mass`, `mass_da`, or `mw`. Example:

```csv
mass
674.2382
965.3336
```

Imports append unique values to the corresponding field. Invalid files leave the previous values unchanged. Files are parsed in the browser; only the resulting numbers are submitted with a search. Limit: 50 KB per file and 20 unique custom masses per type.

## Mass calculations

All calculations are monoisotopic. Proton mass is **1.007276466621 Da** and water mass is **18.010564684 Da**. Peptide neutral mass is the sum of amino acid residue masses plus one terminal water; no water is subtracted between residues.

Choose **Free neutral glycan mass** (default) or **Attached glycan residue mass**. For free glycans, subtract one water from the supplied glycan before adding it to the peptide. For attached residue/composition masses, add the supplied mass directly without another water loss. Library selections convert automatically; custom numbers stay as entered and are interpreted according to the selected type. ExPASy library free masses retain the calculator’s reducing-end convention; library attached increments are explicitly derived as free mass minus the configured water, rather than presented as exact sums of the rounded sugar residue table.

Choose **Positive** (default) or **Negative** ion mode. Charge is always a positive integer magnitude:

- Positive: calculated m/z = neutral mass / z + proton mass; target neutral mass = observed m/z × z − z × proton mass.
- Negative: calculated m/z = neutral mass / z − proton mass; target neutral mass = observed m/z × z + z × proton mass.

Results display absolute Δm/z as the main error. Signed Δm/z = calculated m/z − observed m/z; signed ppm = signed Δm/z / observed m/z × 10⁶. The calculation view separately labels ΔMass as the neutral mass difference in Da. Tolerance can be ppm, m/z, or neutral-mass Da. The closest 25 candidates are ranked by absolute error for each measurement, and all in-tolerance candidates are counted.

The API accepts `glycan_mass_type: free | attached` and `ion_mode: positive | negative`. Legacy `mode` values `user` and `mono_free` map to free glycans; `mono_attached` maps to attached glycans. An explicit glycan mass type takes precedence. Legacy requests now default to positive ions; clients requiring negative ions must specify that mode. The old `corrected_mass` field is removed. `mass` / `neutral_mass` and `target_mass` / `target_neutral_mass` are neutral quantities; `error_da` is signed ΔMass. CSV includes the ion mode, glycan mass type, target neutral mass, absolute/signed Δm/z, signed/absolute ppm, and mass breakdown.

AA residue constants remain formula-derived from NIST isotopes with 11 decimal places and integer prefix sums. Display precision does not imply equivalent experimental accuracy.

This is candidate generation, not validated identification or site localization. It does not model multiple glycans, other modifications, disulfide bonds, adducts other than protons, isotope errors, or MS/MS spectra. N sites must be N, and O sites S/T. Atypical N sites are allowed without a sequon requirement. For other O-linked residue chemistries, the validation requires extension.

## Limits and privacy

Maximum 10,000 residues, 100 measurements, 20 custom glycan masses per type plus library selections. A conservative search-size bound of 1.5 million assignments prevents excessive work; shorten the peptide range or reduce sites if needed. The default maximum length is 60 residues, so longer candidates require increasing it. Inputs are not saved to disk. Only accession import makes an external API request to UniProt. Fonts may load from Google Fonts; system fonts work offline. The local server binds to localhost. For public hosting, the Vercel configuration serves static assets and runs the API as a Python Function; see [DEPLOYMENT.md](DEPLOYMENT.md).

## Verify

```sh
python3 -m unittest discover -s tests -v
```

The Python tests cover reference masses, glycan model conversion, calculation totals, custom masses, isobaric references, real UniProt annotations, missing annotations, unsupported sites, and isoform handling, as well as the original search behavior. `tests/fixtures/P02787-glycosylation.json` is a reduced public UniProt record used for deterministic tests.

### Browser checks

Optional end-to-end checks require Node.js 22+ and Google Chrome. Start the local app as above. In a second terminal, start an isolated headless Chrome session on macOS:

```sh
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
  --headless --disable-gpu --no-first-run --no-default-browser-check \
  --remote-debugging-port=9223 --remote-debugging-address=127.0.0.1 \
  --user-data-dir=/tmp/glycotag-browser-profile about:blank
```

In a third terminal, from the project folder, test the local website:

```sh
node tests/browser_smoke.mjs
```

Or test the public website (this sends sample test requests to the live site):

```sh
GLYCOTAG_URL=https://glycotag.vercel.app node tests/browser_smoke.mjs
```

The checks exercise live UniProt import, library selection, structure dialogs, matching, calculations, custom file import, CSV export, automatic batching, empty-annotation handling, and mobile layout. Screenshots are written to `/tmp/glycotag-*.png`. When finished, press **Ctrl+C** in the headless Chrome terminal; stop the local Python server separately with **Ctrl+C** in its terminal.

Files: `app.py` HTTP server, `engine.py` matching logic and calculation breakdowns, `masses.py` precise amino acid and atomic masses, `uniprot.py` annotation import, `glycans.py` glycan references and masses, `static/` browser interface, `tests/` validation.

Scientific references: [ExPASy PeptideMass](https://web.expasy.org/peptide_mass/peptide-mass-doc.html), [ExPASy masses](https://web.expasy.org/findmod/findmod_masses.html), [UniProt P02787](https://www.uniprot.org/uniprotkb/P02787/entry).

Glycan references: [GlycoMod sugar masses](https://web.expasy.org/glycomod/glycomod_masses.html), [O-GalNAc glycans](https://www.ncbi.nlm.nih.gov/books/NBK579921/), [N-glycans](https://www.ncbi.nlm.nih.gov/books/NBK579964/), [O-GlcNAc](https://www.ncbi.nlm.nih.gov/books/NBK1954/), and [SNFG symbols](https://www.ncbi.nlm.nih.gov/glycans/snfg.html). Drawings are generated from the reference connectivity data in the code; they are not copied images or experimentally inferred structures.

Atomic mass source: [NIST isotope masses](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?isotype=all). Residue formula reference: [University of Washington Proteomics Resource](https://proteomicsresource.washington.edu/resources/data-analysis/masses/). The app derives masses using the NIST values, so final digits can differ from tables based on older atomic constants.

## Public hosting

Vercel deployment files are included: `vercel.json`, `api/index.py`, `.python-version`, `requirements.txt`, and `.vercelignore`. See [DEPLOYMENT.md](DEPLOYMENT.md) for account setup, login, deployment, and public-link verification. Deployment requires authenticating your own Vercel account.
