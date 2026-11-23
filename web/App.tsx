import { useRef, useState } from "react";
import { scoreSample, scoreText, type Report, type Verdict } from "./score";

const samples = [
  {
    id: "helix" as const,
    label: "Clean helix",
    detail: "Ideal backbone, pLDDT 93, low PAE.",
  },
  {
    id: "clash" as const,
    label: "Overconfident clash",
    detail: "Atoms piled together, pLDDT still 92.",
  },
  {
    id: "shift" as const,
    label: "Shifted domains",
    detail: "Local fold kept, domains pulled 16 Å apart.",
  },
];

export function App() {
  const structureText = useRef<string | null>(null);
  const paeText = useRef<string | null>(null);
  const [structureName, setStructureName] = useState<string | null>(null);
  const [paeName, setPaeName] = useState<string | null>(null);
  const [sampleId, setSampleId] = useState<(typeof samples)[number]["id"] | null>(null);
  const [report, setReport] = useState<Report | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  function run(next: () => Report) {
    setBusy(true);
    setError(null);
    window.setTimeout(() => {
      try {
        setReport(next());
      } catch (caught) {
        setReport(null);
        setError(caught instanceof Error ? caught.message : "Could not score that file.");
      } finally {
        setBusy(false);
      }
    }, 30);
  }

  function scoreCurrent(structure: string, name: string) {
    run(() => scoreText(name, structure, paeText.current));
  }

  async function onStructure(file: File) {
    const text = await file.text();
    structureText.current = text;
    setStructureName(file.name);
    setSampleId(null);
    scoreCurrent(text, file.name);
  }

  async function onPae(file: File) {
    const text = await file.text();
    paeText.current = text;
    setPaeName(file.name);
    if (!structureText.current || !structureName) {
      setError("Load a PDB or mmCIF structure before the PAE file.");
      return;
    }
    scoreCurrent(structureText.current, structureName);
  }

  function onSample(id: (typeof samples)[number]["id"]) {
    structureText.current = null;
    paeText.current = null;
    setStructureName(null);
    setPaeName(null);
    setSampleId(id);
    run(() => scoreSample(id));
  }

  return (
    <main className="mx-auto max-w-6xl px-4 py-8 sm:px-8 sm:py-12">
      <header className="max-w-3xl">
        <p className="text-xs font-semibold tracking-[0.22em] text-[#b85a32] uppercase">AlphaFold review</p>
        <h1 className="mt-2 display text-5xl font-medium tracking-tight text-[#1c1712] sm:text-6xl">
          Foldcheck
        </h1>
        <p className="mt-4 max-w-2xl text-lg leading-relaxed text-[#4a4036]">
          A high pLDDT does not mean the fold is right. Foldcheck reads a predicted structure, checks the
          backbone, and estimates how close the coordinates are to a real fold.
        </p>
      </header>

      <div className="mt-10 grid gap-6 lg:grid-cols-[minmax(0,0.9fr)_minmax(0,1.1fr)]">
        <section className="rounded-3xl border border-[#e0d2be] bg-[#f7f1e6] p-5 shadow-[0_20px_50px_rgba(70,42,18,0.05)] sm:p-6">
          <h2 className="display text-2xl">Structure</h2>
          <p className="mt-2 text-sm leading-relaxed text-[#5c5146]">
            Drop an AlphaFold PDB or mmCIF. Files stay in this browser and are not uploaded.
          </p>
          <DropZone
            label="PDB or mmCIF"
            hint={structureName ?? "ranked_0.pdb, model.cif"}
            accept=".pdb,.cif,.mmcif,.ent"
            onFile={(file) => void onStructure(file)}
          />
          <DropZone
            label="Optional PAE or confidence JSON"
            hint={paeName ?? "predicted_aligned_error.json"}
            accept=".json"
            onFile={(file) => void onPae(file)}
          />
          <div className="mt-6">
            <h3 className="text-xs font-semibold tracking-[0.18em] text-[#7a6b5c] uppercase">Try a known case</h3>
            <div className="mt-3 grid gap-2">
              {samples.map((sample) => (
                <button
                  key={sample.id}
                  type="button"
                  onClick={() => onSample(sample.id)}
                  className={`rounded-2xl border px-4 py-3 text-left transition ${
                    sampleId === sample.id
                      ? "border-[#b85a32] bg-[#fff8f3]"
                      : "border-[#e4d7c6] bg-white/50 hover:border-[#c9b59d]"
                  }`}
                >
                  <span className="block font-semibold">{sample.label}</span>
                  <span className="mt-1 block text-sm text-[#5c5146]">{sample.detail}</span>
                </button>
              ))}
            </div>
          </div>
        </section>

        <section className="rounded-3xl border border-[#e0d2be] bg-[#fffaf3] p-5 sm:p-6">
          {busy && <p className="text-[#5c5146]">Reading coordinates…</p>}
          {!busy && error && (
            <div>
              <h2 className="display text-3xl text-[#9d342c]">Could not score</h2>
              <p className="mt-3 leading-relaxed text-[#4a4036]">{error}</p>
            </div>
          )}
          {!busy && !error && !report && (
            <div>
              <h2 className="display text-3xl">Nothing scored yet</h2>
              <p className="mt-3 leading-relaxed text-[#4a4036]">
                Load one of your AlphaFold models, or open a sample. A reliable call needs clean geometry and a
                high, narrow accuracy estimate. A pretty confidence number is not enough.
              </p>
              <dl className="mt-6 grid gap-4 sm:grid-cols-3">
                <Fact title="Geometry" body="Clashes, bond lengths, Ramachandran outliers, and chain breaks." />
                <Fact title="Estimated lDDT" body="How close the coordinates likely are to a true structure." />
                <Fact title="Viability" body="Both of those, with a penalty when pLDDT and geometry disagree." />
              </dl>
            </div>
          )}
          {!busy && !error && report && <Result report={report} />}
        </section>
      </div>
    </main>
  );
}

function Result({ report }: { report: Report }) {
  return (
    <div>
      <div className="flex flex-wrap items-end justify-between gap-3">
        <div>
          <p className="text-xs font-semibold tracking-[0.18em] text-[#7a6b5c] uppercase">{report.name}</p>
          <h2 className={`mt-1 display text-5xl ${verdictClass(report.verdict)}`}>
            {label(report.verdict)}
          </h2>
        </div>
        <p className="text-sm text-[#5c5146]">
          {report.residues} residues · {report.chains} {report.chains === 1 ? "chain" : "chains"} · {report.atoms}{" "}
          atoms
        </p>
      </div>
      <div className="mt-6 grid gap-3 sm:grid-cols-3">
        <Metric label="Viability" value={report.viability} />
        <Metric label="Geometry" value={report.geometryScore} />
        <Metric
          label="Estimated lDDT"
          value={report.predictedLddt}
          detail={`${report.predictedLow.toFixed(2)}–${report.predictedHigh.toFixed(2)}`}
        />
      </div>
      {report.plddt && (
        <div className="mt-6">
          <div className="flex items-baseline justify-between gap-3">
            <h3 className="text-sm font-semibold">pLDDT along the chain</h3>
            <p className="text-sm text-[#5c5146]">{report.plddtMean == null ? "—" : report.plddtMean.toFixed(0)} mean</p>
          </div>
          <div className="mt-2 flex h-8 overflow-hidden rounded-full border border-[#eadfce]">
            {report.plddt.map((value, index) => (
              <span
                key={`${index}-${value.toFixed(1)}`}
                title={`Residue ${index + 1}: ${value.toFixed(0)}`}
                className="h-full flex-1"
                style={{ background: plddtColor(value) }}
              />
            ))}
          </div>
          <p className="mt-2 text-xs text-[#7a6b5c]">Orange is below 50. Blue is 90 or higher.</p>
        </div>
      )}
      <ul className="mt-6 space-y-2">
        {report.reasons.map((reason) => (
          <li key={reason} className="rounded-2xl bg-[#f6efe4] px-4 py-3 text-sm leading-relaxed text-[#2b241d]">
            {reason}
          </li>
        ))}
      </ul>
    </div>
  );
}

function DropZone({
  label,
  hint,
  accept,
  onFile,
}: {
  label: string;
  hint: string;
  accept: string;
  onFile: (file: File) => void;
}) {
  return (
    <label
      className="mt-4 block cursor-pointer rounded-2xl border border-dashed border-[#d3c2ab] bg-white/60 px-4 py-4 hover:border-[#b85a32]"
      onDragOver={(event) => event.preventDefault()}
      onDrop={(event) => {
        event.preventDefault();
        const file = event.dataTransfer.files?.[0];
        if (file) onFile(file);
      }}
    >
      <span className="block font-semibold">{label}</span>
      <span className="mt-1 block text-sm text-[#5c5146]">{hint}</span>
      <input
        className="sr-only"
        type="file"
        accept={accept}
        onChange={(event) => {
          const file = event.target.files?.[0];
          if (file) onFile(file);
        }}
      />
    </label>
  );
}

function Metric({ label, value, detail }: { label: string; value: number; detail?: string }) {
  return (
    <div className="rounded-2xl bg-[#1c1712] px-4 py-3 text-[#f6efe4]">
      <p className="text-xs tracking-[0.14em] text-[#d9c7b2] uppercase">{label}</p>
      <p className="mt-1 display text-3xl">{value.toFixed(2)}</p>
      {detail && <p className="text-xs text-[#d9c7b2]">{detail}</p>}
    </div>
  );
}

function Fact({ title, body }: { title: string; body: string }) {
  return (
    <div className="rounded-2xl bg-[#f6efe4] p-3">
      <dt className="font-semibold">{title}</dt>
      <dd className="mt-1 text-sm leading-relaxed text-[#5c5146]">{body}</dd>
    </div>
  );
}

function label(verdict: Verdict): string {
  if (verdict === "reliable") return "Reliable";
  if (verdict === "reject") return "Reject";
  return "Uncertain";
}

function verdictClass(verdict: Verdict): string {
  if (verdict === "reliable") return "text-[#0e6b52]";
  if (verdict === "reject") return "text-[#9d342c]";
  return "text-[#8d5a12]";
}

function plddtColor(value: number): string {
  if (value >= 90) return "#1d4e89";
  if (value >= 70) return "#3d7ea6";
  if (value >= 50) return "#f0b429";
  return "#e07a3d";
}
