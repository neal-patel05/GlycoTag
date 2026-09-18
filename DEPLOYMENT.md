# Deploy GlycoTag to Vercel

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

Select your personal Hobby scope. Create a new project named `glycotag` (or choose another available name), use the current folder, and retain the settings supplied by `vercel.json`. Framework: Other; output directory: `static`; no build command. No environment variables are required. The CLI uploads directly from this folder; Git integration is optional.

The resulting production alias will be a Vercel-provided HTTPS `.vercel.app` URL. Use the actual alias printed by the successful deployment, not a guessed URL. A generated preview URL may require authentication; verify that the production alias can be opened signed out. If necessary, review Project Settings → Deployment Protection and ensure the production site is public; preview protection can remain enabled.

## Verify before sharing

- Open the production URL in a signed-out / private window.
- Confirm the page, stylesheet and JavaScript load.
- Confirm `/api/glycans` returns JSON and the 21-entry library appears.
- Import `P02787`; verify S51, N432, N491 and N630.
- Search a known peptide and confirm the amino acid precision, proton subtraction, calculation breakdown, and CSV export.
- Verify bad input produces a useful error and that results for one visitor are not shown to another.

## Architecture

`static/` is the only published static directory. `/api/*` requests are rewritten to `api/index.py`, which shares the API logic in `app.py`. Python runs inside a Vercel Function with a 60-second maximum duration. `.python-version` selects Python 3.12. No third-party Python libraries are needed. The original `python3 app.py` local command still works.

`.gitignore` and `.vercelignore` exclude credentials, local Vercel metadata, tests, and development artifacts from publishing. Search POST responses are marked `no-store`. The browser submits large measurement lists in batches of 10 to stay within hosted response limits; exceptionally large direct API responses return a clear 413 error. The app does not intentionally persist submitted sequences or measurements; hosted requests are processed on Vercel and the platform may retain operational metadata. UniProt imports call the public UniProt API.

Search limits remain in place: 10,000 residues, 100 measurements, and a conservative 1.5-million-candidate bound. Vercel plan quotas still apply. If traffic grows, inspect the Vercel Usage page before changing plans or enabling paid features.

## Future changes

Edit and test locally, then run `npx vercel --prod` from the same linked folder. It updates the existing project; the production alias remains stable. Keep `.vercel/` on this computer but do not commit or publish it. The Vercel dashboard shows deployment status, build logs, runtime logs, and rollback controls.

References: [Python API functions](https://vercel.com/docs/functions/runtimes/python/api-directory), [CLI deployments](https://vercel.com/docs/projects/deploy-from-cli), [project configuration](https://vercel.com/docs/project-configuration/vercel-json).
