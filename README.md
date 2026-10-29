# GlyPeptID

**Local app: http://localhost:8000** · **Public address: https://glypeptid.vercel.app**

Source repository: [neal-patel05/GlycoTag](https://github.com/neal-patel05/GlycoTag). The Vercel project is **neals-team/glypeptid**, renamed from glycotag while retaining its project ID. Production deployment uses the Vercel CLI; see [DEPLOYMENT.md](DEPLOYMENT.md).

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

Stopping the local server does **not** stop an existing public deployment.

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

If it is the GlyPeptID Python server you want to stop, note its `PID` and run `kill` with that number. For example, **only if the displayed PID is 12345**:

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

Verify the public address:

```sh
curl --fail --silent --show-error https://glypeptid.vercel.app/api/example
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

Deployments update the existing **glypeptid** project. Local edits do not reach the public website until you deploy. Verify the production alias printed by the CLI after each deployment.

For a preview deployment that does not replace production:

```sh
npx vercel --scope neals-team
```

Open the preview URL printed by the CLI. Preview URLs may require Vercel authentication.

Keep `.vercel/` project linkage locally; do not commit it. The project ID remains unchanged after the rename to `glypeptid`. More details: [DEPLOYMENT.md](DEPLOYMENT.md).

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
4. Set tolerance, peptide length range, and calculation model. Select common glycans, browse the separate **GlyCosmos expanded library**, or supply custom masses. The complete GlyCosmos catalogue includes all glycan types. Select individual records, a page, or all 257,024; choose which protein site types to test separately. Each library has its own clear-selection control.
5. In **Adduct selection**, keep proton-only mode or choose Fiehn Lab adducts. Positive and negative sections each support select all, clear, and individual checkboxes. The active ion mode and measurement charge filter the selected adducts.
6. Run the search. Use **View calculation** for any candidate to see the amino acids, water corrections, glycan contribution, target mass, and signed difference. Download the displayed candidates and calculation totals as CSV.

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

## Expanded GlyCosmos library

The **GlyCosmos complete library** starts collapsed. Click its heading to expand all browsing and selection options. It includes all **257,024 records** in the downloaded catalogue (retrieved 2026-10-26): N-linked, O-linked, other glycans, and unclassified records. Search by accession, IUPAC text, motif, or composition; browse 30 records per page. **Filter by glycan linkage** offers **N-linked only**, **O-linked only**, **Both N-linked and O-linked**, and **All glycans**. **Select all** spans the chosen species and linkage category across every page, regardless of the text search or supported-mass filter. Changing species or linkage clears the current GlyCosmos selection so hidden selections cannot enter a differently filtered search. **Select this page**, individual checkboxes, individual exclusions after select-all, and **Clear selection** are also available. Records without supported masses can still be selected; matching reports selected, supported, and skipped counts.

**Species / organism** filters the catalogue using the source CSV’s **Taxonomy** annotations from [GlyCosmos](https://glycosmos.org/glycans). Type **human / humans** (Homo sapiens, 9606), **bovine / cow / cattle** (Bos taurus, 9913), a scientific name, or an exact taxonomy ID. Suggestions show source names, taxonomy IDs and glycan counts. An exact name or alias selects that organism; partial names match all listed organisms containing the text until you choose a specific suggestion. This is a direct annotation filter, not a taxonomic descendant or inferred-association search. Clear the field to restore all organisms.

The current snapshot has **21,689 glycans with organism annotations**, covering **2,027 taxonomy entries**, including **8,016 human** and **1,755 bovine** glycans. The other **235,335 records have no Taxonomy annotation** and remain available with the species field blank. Reported occurrence is not species exclusivity. With **All glycans**, the species filter includes N-linked, O-linked, other and unclassified structures; combine it with the linkage menu to narrow further. Select-all, individual selections, exclusions, both calculation tools, and CSV provenance retain the species scope. These filters apply to the GlyCosmos library; common-library and custom masses remain separate selections.

Cards show the reported organisms and source IUPAC name (or a clearly labelled accession when none is supplied), a structure image where supplied (local archive or official hosted image), and residue abbreviations with full subunit names. **Structure & subunit names** opens a larger view, mass details, source links and WURCS notation. Sugar labels come from 8,961 single-residue GlyCosmos records plus the GlyCosmos residue-name table; unresolved residue chemistry remains explicitly unresolved. Repeating structures label counts as repeat-unit counts. Motif annotations are not represented as unique glycan names.

The official image archive supplies **192,581 structure images** for this catalogue. They are installed in the local image cache and served without external requests. The remaining **64,443 records** have no image in that archive and are clearly marked. No structure or missing name is invented.

Linkage filters use the catalogue’s **Motif name(s)** accessions, mapped to N/O categories using [GlyCosmos GlycoNAVI motifs](https://glycosmos.org/glyconavi/motifs). The current snapshot contains **22,667 N-linked**, **8,850 O-linked**, and **31,517 combined** annotated records. These are motif-based categories, not experimentally verified protein attachments. The mapping includes N-glycan core/high-mannose/hybrid/complex motifs; O-glycan cores 1–7 (including fuzzy forms), T/Tn antigens, sialylated T/Tn, O-fucose and O-mannose variants. The Ser/Thr categories were checked against the [GlyCosmos motif table](https://beta.glycosmos.org/glyconavi/motifs?page=7&q%5Bs%5D=aglycon+asc); sialyl-Tn is also supported by [GlyCosmos glycosylation records](https://glycosmos.org/glycans/G36123IU?taxonomy_page=2). No linkage is inferred from mass or a generic sugar name such as N-acetyllactosamine. Unannotated records remain available under **All glycans**. **Match glycans to compatible protein sites** can restrict matching to N sites, O sites, or both; it cannot override a glycan’s linkage. The peptide engine assigns GlyCosmos N-linked records only to Asn sites and O-linked records only to Ser/Thr sites. Glycans without an N/O source annotation are excluded from peptide-site matching. They remain available for browsing and standalone glycan/adduct matching, which does not assign a protein attachment site. Custom masses use the N/O field where they were entered. Sites with no compatible selected masses are skipped and listed in the results; if no compatible sites remain, the search reports an explanatory error.

A mass index and bounded nearest-candidate search support selecting the full library without materializing the peptide × glycan product. Isobaric records share calculation work. Results report the total number of references sharing a mass and display up to 20 accessions with their source structure images. The first reference is shown in each result row; expand **Show more reference structures** for other references sharing the mass. Calculation details also show the images. **Structure & subunit names** opens the full record even when that glycan has never been browsed. Missing images are explicitly marked, and mass matches do not uniquely identify a structure. CSV includes the reference count, linkage filter, and compact full-library selection/exclusions.

**116,089 records** currently have supported ExPASy composition masses. Of these, **28,942** have an N/O motif annotation and are eligible for peptide-site matching when compatible sites are enabled; the remaining **87,147** have no N/O annotation. Results and CSV coverage distinguish unavailable masses (`skipped`) from supported records excluded by linkage (`linkage_skipped`, split into `unclassified` and `incompatible`). `matched` counts the eligible records. Each site type has its own mass/reference groups, so even identical-mass N/O isomers cannot appear under the wrong attachment site. Every record remains browsable; unsupported chemistry, reduced ends, repeats, cyclic structures, incomplete compositions, and source-mass disagreements are explicitly unavailable for matching. We do not substitute GlyCosmos's two-decimal mass as a precision ExPASy result. A conservative WURCS parser maps fixed supported compositions to the ten ExPASy inputs. It also checks calculated masses against the source's rounded mass when available. This is a restricted composition conversion, not a general WURCS chemistry engine.

The ten ExPASy underivatized monoisotopic residue values and reducing-end convention were checked directly against GlycanMass on 2026-10-26. Results are saved in `tests/fixtures/expasy-expanded.json`. Masses are calculated locally with those constants; no bulk requests are sent to ExPASy. Free/attached conversion is the same as the common library.

Source: [GlyCosmos full catalogue](https://download.glycosmos.org/glycosmos_glycans_list.csv), [GlyCosmos downloads](https://glycosmos.org/download), licensed CC BY 4.0. `data/glycosmos.metadata.json` records counts, provenance, source SHA-256 and coverage exclusions. The compressed SQLite snapshot is bundled; first use extracts a read-only database to the OS temporary directory. No network is needed to browse or search it. The uncompressed build database is ignored by Git and deployments.

To rebuild locally from a newly downloaded full CSV (then restart the server):

```sh
python3 scripts/build_glycosmos.py /path/to/glycosmos_glycans_list.csv
```

When rebuilding an existing download, add `--retrieved YYYY-MM-DD` with its original download date. The builder retains and indexes organism names and taxonomy IDs alongside the glycan records. `GET /api/glycosmos?species=human&linkage=both` returns filtered cards, matching organisms, suggestions, and species-scoped linkage counts. Analysis requests carry the same filter in `glycosmos_selection.species`.

The local images can be rebuilt from the public archive:

```sh
curl -L --fail https://image.glycosmos.org/snfg-png.zip -o /tmp/snfg-png.zip
python3 scripts/build_glycosmos_images.py /tmp/snfg-png.zip
```

The download is about 468 MB; the local indexed image cache is about 624 MB and is ignored by Git/deployment. Its source hash and count are in `data/glycosmos-images.metadata.json`. The small `glycosmos-image-ids.json.gz` manifest is included in Git and deployments. When the local archive is absent, image requests redirect to GlyCosmos’s official HTTPS image service; the catalogue, species annotations, and mass calculations remain bundled. Hosted image availability depends on that external service. Restart the app after replacing either catalogue or images. HTML/JS/API responses use `Cache-Control: no-store` to prevent stale interface assets. If the library endpoint returns “Not found,” restart the existing Python server; refreshing HTML alone cannot update an old Python process.

## Adduct library

All **31 positive** and **15 negative** entries from the [Fiehn Lab calculator](https://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator/) are included, using the numeric column of its corrected 2020 spreadsheet. This corrects the web table's M+H+2Na offset and 3M-H sign. Other entries preserve the published precision and labels, including the source's unusual IsoProp+Na+H entry; they are tabulated hypotheses, not newly validated chemical species. Pure proton ions retain the app's more precise proton constant.

Both **Positive ion mode (31)** and **Negative ion mode (15)** sections are visible and expanded on initial load, independently of the proton-only setting. Select individual entries or all entries in either section; selecting an entry or select-all activates library mode. Selections persist across polarity changes; only the active polarity is searched. Fixed-charge entries are evaluated only for matching measurement charges. The default proton-only option supports arbitrary charges up to 100. If no selected adduct is compatible with a measurement, it gets a clear no-compatible-adduct result instead of a fabricated match.

For each adduct, **m/z = (n × M + shift) / z** and **target M = (observed m/z × z − shift) / n**, where M is the neutral mass of one glycopeptide, n is the multimer count, and z is charge magnitude. Dimer/trimer neutral-mass tolerances refer to one monomer. Rankings compare absolute m/z errors across all compatible selected adducts. Result rows, calculation details, and CSV include the adduct, shift, multimer count, and the candidate-specific target. A result with several adduct hypotheses has no single neutral target.

### One best glycan signal assignment

Click **Identify best glycan signal** in Adduct selection to compare **glycan-only** against **adduct-associated glycan** signals using the current measurements, glycan selections, polarity, adduct selections and tolerance. This action requires no protein sequence or sites. It reports exactly one assignment per unique **m/z and charge** in both the interface and its CSV export.

1. Calculate **MW = (m/z × z) − z**, retaining the requested literal formula in either ion mode.
2. Compare all selected **free neutral glycan masses** under the two explanations. Glycan-only means monomeric protonation/deprotonation with no additional adduct. This baseline always competes, including when no selected additional adduct supports the charge. Equivalent proton-only library entries are deduplicated. Library masses use their free values; custom attached masses are converted by adding water.
3. Choose the smallest absolute m/z error among explanations within tolerance. If none fit, report the closest explanation and label it **outside tolerance**, rather than claiming an identification. Numerical ties within 1e-9 m/z are flagged; one deterministic representative is shown, preferring glycan-only, then lower glycan mass and adduct ID. Isobaric glycan reference structures remain ambiguous and are grouped under that single assignment.
4. Show the winning category, glycan IDs and structures, **difference = assigned glycan MW − calculated MW**, predicted m/z, errors and tolerance status. An additional adduct is populated only for an adduct-associated winner. The CSV has one row per measurement, with `assignment_type`, ion form, optional adduct, tie count and calculation values.

Subtracting z assumes 1 Da per charge, so calculated MW is labelled an approximation; it is not a universal exact neutral-mass conversion. Ion comparisons use the precise proton mass and Fiehn offsets. The expected difference for glycan mass G is **z − shift − (n − 1) × G**, derived from **z × m/z = n × G + shift**. Difference error divided by z gives signed m/z error; neutral-mass Da tolerance uses difference error divided by n. For example, m/z 1022.989218 at charge 1 gives calculated MW 1021.989218; for G = 1000, the difference is −21.989218 Da, agreeing with `[M+Na]+` (shift +22.989218 Da).

The best assignment searches across all selected glycan masses and compatible ion models, using binary-search neighbors for each model. It does not restrict adduct testing to the glycan nearest the approximate MW. This allows an exact sodium match at G = 1000 to beat a closer approximate mass at G = 1022. These are mutually exclusive reporting categories, not a claim that glycans cannot carry adducts chemically.

The dedicated API is `POST /api/glycan-adducts`; it accepts the same measurement, glycan and adduct fields as the peptide search, but does not use protein sequence, sites, peptide lengths or `glycosmos_apply_to`. Each result returns one `best_assignment` and a `tied_assignment_count`; the former multi-candidate/adduct lists are no longer returned. Species and N/O GlyCosmos catalogue filters and exclusions still apply. Unsupported glycan masses are skipped and reported.

## Mass calculations

All calculations are monoisotopic. Proton mass is **1.007276466621 Da** and water mass is **18.010564684 Da**. Peptide neutral mass is the sum of amino acid residue masses plus one terminal water; no water is subtracted between residues.

Choose **Free neutral glycan mass** (default) or **Attached glycan residue mass**. For free glycans, subtract one water from the supplied glycan before adding it to the peptide. For attached residue/composition masses, add the supplied mass directly without another water loss. Library selections convert automatically; custom numbers stay as entered and are interpreted according to the selected type. ExPASy library free masses retain the calculator’s reducing-end convention; library attached increments are explicitly derived as free mass minus the configured water, rather than presented as exact sums of the rounded sugar residue table.

Choose **Positive** (default) or **Negative** ion mode. Charge is always a positive integer magnitude:

In proton-only mode:

- Positive: calculated m/z = neutral mass / z + proton mass; target neutral mass = observed m/z × z − z × proton mass.
- Negative: calculated m/z = neutral mass / z − proton mass; target neutral mass = observed m/z × z + z × proton mass.

Results display absolute Δm/z as the main error. Signed Δm/z = calculated m/z − observed m/z; signed ppm = signed Δm/z / observed m/z × 10⁶. The calculation view separately labels ΔMass as the neutral mass difference in Da. Tolerance can be ppm, m/z, or neutral-mass Da. The closest 25 candidates are ranked by absolute error for each measurement, and all in-tolerance candidates are counted.

The API accepts `glycan_mass_type: free | attached` and `ion_mode: positive | negative`. Legacy `mode` values `user` and `mono_free` map to free glycans; `mono_attached` maps to attached glycans. An explicit glycan mass type takes precedence. Legacy requests now default to positive ions; clients requiring negative ions must specify that mode. The old `corrected_mass` field is removed. `mass` / `neutral_mass` and `target_mass` / `target_neutral_mass` are neutral quantities; `error_da` is signed ΔMass. CSV includes the ion mode, glycan mass type, target neutral mass, absolute/signed Δm/z, signed/absolute ppm, and mass breakdown.

`GET /api/adducts` returns valid adduct IDs. Supply `selected_adducts: ["pos_15"]` to test sodium; omit the field for arbitrary-charge proton-only behavior. An explicit empty selection is rejected. `GET /api/glycosmos?q=HexNAc&page=1&supported=1` browses the expanded library. Supply `glycosmos_selection: {"ids":["G00039MO"]}` for individual records, or `glycosmos_selection: {"all":true,"excluded":[]}` for the full catalogue. `glycosmos_apply_to: "both" | "N" | "O"` restricts compatible protein site types without overriding glycan linkage. Legacy `glycosmos_glycans` requests remain supported and reject incompatible or unannotated protein-site assignments. `GET /api/glycosmos/image/G00039MO` serves the local PNG. Client-supplied library masses are never trusted. Each result row has its own `target_neutral_mass`; the measurement-level target is `null` when there are multiple or no compatible adduct hypotheses.

AA residue constants remain formula-derived from NIST isotopes with 11 decimal places and integer prefix sums. Display precision does not imply equivalent experimental accuracy.

This is candidate generation, not validated identification or site localization. It does not model multiple glycans, other modifications, disulfide bonds, isotope errors, or MS/MS spectra. N sites must be N, and O sites S/T. Atypical N sites are allowed without a sequon requirement. For other O-linked residue chemistries, the validation requires extension.

## Limits and privacy

Maximum 10,000 residues, 100 measurements, 20 custom glycan masses per type plus library selections, including the entire GlyCosmos catalogue. Work is bounded at 1.5 million peptide/site combinations and 1.5 million peptide × compatible-adduct × measurement queries. Glycan mass alternatives use binary searches rather than expanding that product; shorten the peptide range, reduce sites, adducts or measurements if needed. The default maximum length is 60 residues, so longer candidates require increasing it. Inputs are not saved to disk. Accession imports contact UniProt, and hosted structure images contact GlyCosmos; the GlyCosmos and adduct catalogues are bundled locally. Fonts may load from Google Fonts; system fonts work offline. The local server binds to localhost. For public hosting, the Vercel configuration serves static assets and runs the API as a Python Function; see [DEPLOYMENT.md](DEPLOYMENT.md).

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
  --user-data-dir=/tmp/glypeptid-browser-profile about:blank
```

In a third terminal, from the project folder, test the local website:

```sh
node tests/browser_smoke.mjs
node tests/browser_expanded.mjs
node tests/browser_glycan_adducts.mjs
node tests/browser_species.mjs
node tests/browser_linkage.mjs
```

This command can test the public website (it sends sample requests to that site):

```sh
GLYPEPTID_URL=https://glypeptid.vercel.app node tests/browser_smoke.mjs
```

The checks exercise live UniProt import, library selection, structure dialogs, matching, calculations, custom file import, CSV export, automatic batching, empty-annotation handling, and mobile layout. The expanded checks cover the complete GlyCosmos catalogue, unavailable masses, all adduct controls, sodium/trimer matching, charge filtering, and CSV provenance. Each script opens and closes its own temporary browser tab. Set `GLYPEPTID_URL` to test another local port; `CHROME_DEBUG_URL` selects another debugger address. Screenshots are written to `/tmp/glypeptid-*.png`. When finished, press **Ctrl+C** in the headless Chrome terminal; stop the local Python server separately with **Ctrl+C** in its terminal.

Files: `app.py` HTTP server, `engine.py` peptide matching and calculation breakdowns, `glycan_adducts.py` single glycan signal assignment, `masses.py` precise amino acid and atomic masses, `uniprot.py` annotation import, `glycans.py` common glycan references, `glycosmos.py` expanded catalogue, `adducts.py` ion hypotheses, `data/` bundled scientific data, `static/` browser interface, `tests/` validation.

Scientific references: [ExPASy PeptideMass](https://web.expasy.org/peptide_mass/peptide-mass-doc.html), [ExPASy masses](https://web.expasy.org/findmod/findmod_masses.html), [UniProt P02787](https://www.uniprot.org/uniprotkb/P02787/entry).

Glycan references: [GlycoMod sugar masses](https://web.expasy.org/glycomod/glycomod_masses.html), [O-GalNAc glycans](https://www.ncbi.nlm.nih.gov/books/NBK579921/), [N-glycans](https://www.ncbi.nlm.nih.gov/books/NBK579964/), [O-GlcNAc](https://www.ncbi.nlm.nih.gov/books/NBK1954/), and [SNFG symbols](https://www.ncbi.nlm.nih.gov/glycans/snfg.html). Drawings are generated from the reference connectivity data in the code; they are not copied images or experimentally inferred structures.

Atomic mass source: [NIST isotope masses](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?isotype=all). Residue formula reference: [University of Washington Proteomics Resource](https://proteomicsresource.washington.edu/resources/data-analysis/masses/). The app derives masses using the NIST values, so final digits can differ from tables based on older atomic constants.

## Public hosting

Vercel deployment files are included: `vercel.json`, `api/index.py`, `.python-version`, `requirements.txt`, and `.vercelignore`. See [DEPLOYMENT.md](DEPLOYMENT.md) for account setup, login, deployment, and public-link verification. Deployment requires authenticating your own Vercel account.
