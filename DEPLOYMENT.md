# GlyPeptID hosting

Production address: **https://glypeptid.vercel.app**. Source: [neal-patel05/GlycoTag](https://github.com/neal-patel05/GlycoTag), branch `main`. The existing Vercel project was renamed from `glycotag` to `glypeptid` on 2026-10-29, preserving project ID `prj_1figRQXV9tAAlJHbTm23ZG30TsMB` and team `neals-team`. Release commit label: **new url name + adducts addition + complete glycan library**.

The release includes all 257,024 catalogue records, species and N/O filters, all 46 adducts, structure references, and enforced N/O peptide-site compatibility. The 21 MB compressed catalogue and a small image-ID manifest are committed and deployed. The optional 624 MB local image archive is excluded: hosted image requests redirect to GlyCosmos's official image service. No paid storage service is required. The pre-deployment Python suite has 83 passing tests.

## Historical deployment record

Published and verified on September 17, 2026.

- Public URL: https://glycotag.vercel.app
- Project dashboard: https://vercel.com/neals-team/glycotag
- Vercel scope / project: `neals-team/glycotag`
- Verified production deployment: `dpl_91Uxx4XaM9LaQZu3CcC3M8T6hFCC`
- CLI used: Vercel 59.20.0; Python runtime: 3.12

The production alias returned HTTP 200 without authentication. Live browser checks passed for UniProt import, glycan structures, precise mass calculations, target/proton handling, CSV export, 21-measurement automatic batching, input validation, and mobile layout. All 29 Python tests passed.

## Account and authorization

1. Go to https://vercel.com/signup and sign up with your preferred provider or email. Complete any verification directly with Vercel.
2. For a personal, non-commercial project, choose Hobby. Do not start a paid plan or Pro trial unless you intend to use it. See https://vercel.com/docs/plans/hobby for eligibility and quotas.
3. No GitHub repository, purchased domain, payment details, or database is needed for this app's initial Hobby deployment.
4. Authenticate the official CLI in the terminal. Complete the browser login yourself; do not paste passwords, session cookies, or tokens into chat.

```sh
cd "/Users/nealpatel/neal paper"
npx vercel login
```

## Publish

```sh
npx vercel --prod
```

For this existing setup, use scope `neals-team`, retain the linked project ID, and use the settings supplied by `vercel.json`. Framework override: Other; output directory: `static`; no build command. No environment variables are required. The CLI uploads directly from this folder; Git integration is optional.

The resulting production alias will be a Vercel-provided HTTPS `.vercel.app` URL. Use the actual alias printed by the successful deployment, not a guessed URL. A generated preview URL may require authentication; verify that the production alias can be opened signed out. If necessary, review Project Settings → Deployment Protection and ensure the production site is public; preview protection can remain enabled.

## Verify before sharing

- Open the production URL in a signed-out / private window.
- Confirm the page, stylesheet and JavaScript load.
- Confirm `/api/glycans` returns JSON and the 21-entry library appears.
- Confirm `/api/glycosmos` returns 257,024 records; species `human` returns 8,016 and `bovine` returns 1,755 before linkage filtering.
- Confirm `/api/adducts` lists 31 positive and 15 negative ions.
- Open `/api/glycosmos/image/G00039MO` and verify the structure image loads through the official image service.
- Verify N-linked glycans cannot be assigned to O sites, and vice versa.
- Import `P02787`; verify S51, N432, N491 and N630.
- Search a known peptide and confirm the amino acid precision, proton subtraction, calculation breakdown, and CSV export.
- Verify bad input produces a useful error and that results for one visitor are not shown to another.

## Architecture

`static/` is the only published static directory. `/api/*` requests are rewritten to `api/index.py`, which shares the API logic in `app.py`. Python runs inside a Vercel Function with a 60-second maximum duration. `.python-version` selects Python 3.12. No third-party Python libraries are needed. The original `python3 app.py` local command still works.

`.gitignore` and `.vercelignore` exclude credentials, local Vercel metadata, tests, and development artifacts from publishing. Search POST responses are marked `no-store`. The browser submits large measurement lists in batches of 10 to stay within hosted response limits; exceptionally large direct API responses return a clear 413 error. The app does not intentionally persist submitted sequences or measurements; hosted requests are processed on Vercel and the platform may retain operational metadata. UniProt imports call the public UniProt API.

Search limits remain in place: 10,000 residues, 100 measurements, and a conservative 1.5-million-candidate bound. Vercel plan quotas still apply. If traffic grows, inspect the Vercel Usage page before changing plans or enabling paid features.

## Future changes

Edit and test locally, push the reviewed changes to GitHub, then run `npx vercel --prod --scope neals-team` from the same linked folder. Verify `glypeptid.vercel.app` points to the new deployment. Keep `.vercel/` on this computer but do not commit or publish it. The Vercel dashboard shows deployment status, build logs, runtime logs, and rollback controls.

References: [Python API functions](https://vercel.com/docs/functions/runtimes/python/api-directory), [CLI deployments](https://vercel.com/docs/projects/deploy-from-cli), [project configuration](https://vercel.com/docs/project-configuration/vercel-json).
