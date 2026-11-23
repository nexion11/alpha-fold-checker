# Foldcheck

Foldcheck reviews an AlphaFold model. A high pLDDT can still sit on a wrong fold, so the app reports geometry and estimated accuracy separately.

Load a PDB or mmCIF, and optionally an AlphaFold PAE JSON. The score stays in your browser.

- **Reliable** means the backbone is clean and the accuracy estimate is high and narrow.
- **Uncertain** means the signals disagree, the PAE says domains may be misplaced, or confidence is missing.
- **Reject** means serious geometry problems or a low accuracy estimate.

Estimated lDDT is a reading of geometry, pLDDT, and PAE. It is not a comparison to a crystal structure. A wrong fold with high pLDDT, low PAE, and clean local geometry looks the same as a correct one.

## Run locally

```bash
npm install
npm run dev
```

Open the local address Vite prints. `npm run check` scores the built-in helix, clash, and domain-shift cases.

## What the samples show

- **Clean helix** — ideal geometry, pLDDT 93, low PAE.
- **Overconfident clash** — atoms forced together while pLDDT stays at 92.
- **Shifted domains** — the chain is pulled apart and the long-range PAE jumps.

`src/foldcheck` has matching Python helpers for parsing and geometry. The page is the scorer.
