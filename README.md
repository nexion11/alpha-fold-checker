# Alpha Fold Checker

Alpha Fold Checker reviews an AlphaFold prediction. A high pLDDT can still sit on a wrong fold, so the app reports physical geometry and estimated accuracy separately.

Drop a PDB or mmCIF, and optionally an AlphaFold PAE JSON. The file stays in your browser.

## Run

```bash
npm install
npm run dev
```

Open the address Vite prints. `npm run check` scores the three built-in cases.

## Verdicts

- **Reliable** — the backbone is clean and the accuracy estimate is high and narrow.
- **Uncertain** — the signals disagree, the estimate is wide, or long-range PAE says the domains may be misplaced.
- **Reject** — the geometry is badly broken, or the accuracy estimate is low.

## Samples

- **Clean helix** — ideal geometry, pLDDT 93, low PAE.
- **Overconfident clash** — atoms forced together while pLDDT stays at 92.
- **Shifted domains** — the chain is pulled apart and the long-range PAE jumps.

Estimated lDDT is a reading of geometry, pLDDT, and PAE. It is not a comparison to an experimental structure. A wrong fold with high pLDDT, low PAE, and clean local geometry looks the same as a correct one.
