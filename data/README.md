# Scientific data

`glycosmos.sqlite3.gz` is a local snapshot of all 257,024 rows in the
[GlyCosmos glycan export](https://download.glycosmos.org/glycosmos_glycans_list.csv),
downloaded 2026-10-26. Attribution: **GlyCosmos Portal**, licensed
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The snapshot preserves accession, IUPAC, WURCS, motifs, taxonomy annotations and reported mass.
Calculated composition/mass fields and coverage notes are additions by GlyPeptID.
`glycosmos.metadata.json` contains the source checksum and coverage counts.
The source export was listed as updated October 19, 2026 on the download page.

The 116,089 supported composition masses use the ten underivatized,
monoisotopic inputs of [ExPASy GlycanMass](https://web.expasy.org/glycanmass/).
These are locally calculated masses using verified ExPASy constants, not
257,024 individually submitted calculator results. Unsupported records
remain browsable and are excluded from matching. See `glycosmos.py` for
the deliberately restricted conversion and `tests/fixtures/expasy-expanded.json`
for the independent live calculator check.

`adducts.json` records all 46 entries from the **Adducts** worksheet of
[Fiehn Lab's corrected 2020 spreadsheet](https://fiehnlab.ucdavis.edu/images/files/software/ESI-MS-adducts-2020.xls).
It retains source labels, charges, multimer counts and numeric offsets,
including the corrected M+H+2Na and 3M-H values. `adducts.py` applies the
engine's full-precision proton constant to pure proton ions. This snapshot
does not include the separate NIST frequency-list worksheet, which is not
the calculator's positive/negative adduct table.

Rebuild the GlyCosmos snapshot with `python3 scripts/build_glycosmos.py SOURCE.csv`.
No external library is needed to build or read it. The 21 MB compressed
snapshot is the portable artifact; the 128 MB uncompressed build file is
ignored. At runtime it is extracted once into the OS temporary directory
and opened read-only. Restart the app after replacing a snapshot.


The local `glycosmos-images.sqlite3` cache contains all **192,581** images
from [GlyCosmos's SNFG PNG archive](https://image.glycosmos.org/snfg-png.zip),
with source SHA-256 in `glycosmos-images.metadata.json`. The source archive
contains no image for the remaining 64,443 records. Rebuild the cache with
`python3 scripts/build_glycosmos_images.py /path/to/snfg-png.zip`.
This 624 MB cache is local only and ignored by Git/deployment. The small
`glycosmos-image-ids.json.gz` manifest preserves image availability in hosted
installs. Without the local cache, the API redirects known image IDs to
`https://image.glycosmos.org/snfg/png/ACCESSION`, the official image service
listed in the [GlyCosmos help](https://glycosmos.org/help). Missing image IDs
return 404. Hosted image display requires that external service.

`residue-units.json` maps 8,961 single-residue WURCS definitions to their
source IUPAC labels and accession references. `build_glycosmos.py` rebuilds
it from the catalogue. `residue-names.json` transcribes the official
[GlyCosmos residue-name table](https://pages.glycosmos.org/composition2wurcs/public/residue_list.html).
These files provide display labels independently of ExPASy mass support;
unresolved names and variable repeat-unit counts are explicitly labelled.
