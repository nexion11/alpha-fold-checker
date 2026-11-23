export type Vec3 = [number, number, number];

export type Atom = {
  name: string;
  element: string;
  resname: string;
  chain: string;
  resseq: number;
  icode: string;
  xyz: Vec3;
  bfactor: number;
};

export type Residue = {
  resname: string;
  chain: string;
  resseq: number;
  icode: string;
  atoms: Record<string, Atom>;
};

export type Verdict = "reliable" | "uncertain" | "reject";

export type Report = {
  name: string;
  residues: number;
  chains: number;
  atoms: number;
  verdict: Verdict;
  viability: number;
  geometryScore: number;
  predictedLddt: number;
  predictedLow: number;
  predictedHigh: number;
  plddtMean: number | null;
  plddt: number[] | null;
  paeMean: number | null;
  paeLocal: number | null;
  paeLong: number | null;
  clashes: number;
  bondOutlierFraction: number | null;
  angleOutlierFraction: number | null;
  ramaOutlierFraction: number | null;
  chainBreaks: number;
  reasons: string[];
};

const BOND_IDEAL: Record<string, number> = {
  "N|CA": 1.458,
  "CA|C": 1.525,
  "C|O": 1.231,
  "CA|CB": 1.53,
};
const PEPTIDE_BOND = 1.329;
const BOND_TOLERANCE = 0.12;
const ANGLE_IDEAL: Record<string, number> = {
  "N|CA|C": 111.2,
  "CA|C|O": 120.8,
};
const PEPTIDE_ANGLE_CA_C_N = 116.2;
const PEPTIDE_ANGLE_C_N_CA = 121.7;
const ANGLE_TOLERANCE = 10;
const CA_DISTANCE = 3.8;
const HEAVY_CLASH = 2.2;
const CA_ONLY_CLASH = 3;
const CHAIN_BREAK_MIN = 2.9;
const CHAIN_BREAK_MAX = 4.5;
const MAX_ATOMS = 80_000;
const MAX_RESIDUES = 3_000;

const HELIX: [number, number] = [-57, -47];
const CB_DIHEDRAL = 122.6;
const CB_ANGLE = 110.4;

export function scoreText(name: string, structureText: string, paeText: string | null): Report {
  const atoms = parseStructure(name, structureText);
  const residues = groupResidues(atoms);
  if (residues.length < 2) {
    throw new Error("This file has fewer than two residues, so there is no fold to score.");
  }
  if (residues.length > MAX_RESIDUES) {
    throw new Error(`This structure has ${residues.length} residues. Foldcheck scores up to ${MAX_RESIDUES}.`);
  }
  const confidence = paeText ? parseConfidence(paeText) : { plddt: null, pae: null };
  return scoreResidues(name, atoms.length, residues, confidence.plddt, confidence.pae);
}

export function scoreSample(kind: "helix" | "clash" | "shift"): Report {
  const length = 40;
  const built = buildPeptide(length, HELIX[0], HELIX[1]);
  if (kind === "helix") {
    assignBfactor(built, 93);
    const pae = blockPae(length, 3.2, 3.6);
    return scoreResidues("clean helix", built.length, groupResidues(built), null, pae);
  }
  if (kind === "clash") {
    const residues = groupResidues(built);
    const target = residues[3].atoms.CA.xyz;
    const victim = residues[15];
    for (const atom of Object.values(victim.atoms)) {
      atom.xyz = [target[0] + 0.15, target[1] + 0.1, target[2] - 0.05];
    }
    assignBfactor(built, 92);
    return scoreResidues("overconfident clash", built.length, groupResidues(built), null, null);
  }
  const residues = groupResidues(built);
  for (const residue of residues) {
    if (residue.resseq > length / 2) {
      for (const atom of Object.values(residue.atoms)) {
        atom.xyz = [atom.xyz[0] + 16, atom.xyz[1], atom.xyz[2]];
      }
    }
  }
  assignBfactor(built, 88);
  const pae = blockPae(length, 3.4, 22);
  return scoreResidues("shifted domains", built.length, groupResidues(built), null, pae);
}

function scoreResidues(
  name: string,
  atomCount: number,
  residues: Residue[],
  plddtOverride: number[] | null,
  pae: number[][] | null,
): Report {
  const geometry = analyzeGeometry(residues);
  const plddt = resolvePlddt(residues, plddtOverride);
  let paeSummary = summarizePae(pae);
  const reasons: string[] = [];
  if (pae && pae.length !== residues.length) {
    reasons.push(
      `The PAE matrix has ${pae.length} rows and the structure has ${residues.length} residues, so the PAE was left unused.`,
    );
    paeSummary = emptyPae();
    pae = null;
  }

  const plddtMean = plddt ? average(plddt) : null;
  let predicted = estimateAccuracy(geometry, plddtMean, paeSummary);
  const geometryScore = geometryScoreOf(geometry);
  let viability = clamp(0.5 * geometryScore + 0.5 * predicted.point, 0, 1);
  const overconfident = plddtMean != null && plddtMean >= 85 && geometryScore < 0.6;
  if (overconfident) viability = clamp(viability - 0.1, 0, 1);

  const spread = predicted.high - predicted.low;
  let verdict: Verdict;
  if (geometryScore < 0.45 || predicted.point < 0.5 || viability < 0.45) verdict = "reject";
  else if (
    geometryScore >= 0.85 &&
    predicted.point >= 0.8 &&
    spread <= 0.25 &&
    (plddtMean == null || plddtMean >= 70)
  ) {
    verdict = "reliable";
  } else verdict = "uncertain";

  if (plddtMean == null && paeSummary.mean == null && verdict === "reliable") {
    verdict = "uncertain";
    reasons.push("No pLDDT or PAE was found, so the accuracy estimate uses geometry alone.");
  }

  pushGeometryReasons(reasons, geometry, residues.length);
  if (plddtMean != null) {
    reasons.push(
      plddtMean >= 70
        ? `Mean pLDDT is ${plddtMean.toFixed(0)}, which AlphaFold treats as a confident local fold.`
        : `Mean pLDDT is ${plddtMean.toFixed(0)}, below the usual 70 line for a confident fold.`,
    );
  } else {
    reasons.push("No pLDDT was found in the B-factor column or a confidence file.");
  }
  if (predicted.domainGap) {
    reasons.push(
      `Long-range PAE is ${paeSummary.longMean?.toFixed(0)} Å against ${paeSummary.localMean?.toFixed(0)} Å locally. Each piece may be folded, but the arrangement between them is uncertain.`,
    );
  } else if (pae == null) {
    reasons.push("No PAE matrix was used. A wrong domain arrangement can still look clean.");
  }
  if (overconfident) {
    reasons.push("Backbone geometry is poor while pLDDT stays high, so the confidence looks overconfident.");
  }
  reasons.push(
    "This estimates physical plausibility and coordinate accuracy. It does not say whether the protein works in a cell.",
  );

  const chains = new Set(residues.map((residue) => residue.chain)).size;
  return {
    name,
    residues: residues.length,
    chains,
    atoms: atomCount,
    verdict,
    viability,
    geometryScore,
    predictedLddt: predicted.point,
    predictedLow: predicted.low,
    predictedHigh: predicted.high,
    plddtMean,
    plddt,
    paeMean: paeSummary.mean,
    paeLocal: paeSummary.localMean,
    paeLong: paeSummary.longMean,
    clashes: geometry.clashes,
    bondOutlierFraction: geometry.bondOutlierFraction,
    angleOutlierFraction: geometry.angleOutlierFraction,
    ramaOutlierFraction: geometry.ramaOutlierFraction,
    chainBreaks: geometry.chainBreaks,
    reasons,
  };
}

type Geometry = {
  clashes: number;
  bondOutlierFraction: number | null;
  angleOutlierFraction: number | null;
  ramaOutlierFraction: number | null;
  chainBreaks: number;
  residueCount: number;
};

function analyzeGeometry(residues: Residue[]): Geometry {
  return {
    clashes: countClashes(residues),
    bondOutlierFraction: fraction(bondChecks(residues)),
    angleOutlierFraction: fraction(angleChecks(residues)),
    ramaOutlierFraction: fraction(ramaChecks(residues)),
    chainBreaks: caTrace(residues).breaks,
    residueCount: residues.length,
  };
}

function geometryScoreOf(geometry: Geometry): number {
  const penalty =
    Math.min(geometry.clashes / Math.max(geometry.residueCount, 1) / 0.15, 1) * 0.5 +
    (geometry.bondOutlierFraction == null ? 0 : Math.min(geometry.bondOutlierFraction / 0.1, 1) * 0.2) +
    (geometry.angleOutlierFraction == null ? 0 : Math.min(geometry.angleOutlierFraction / 0.15, 1) * 0.1) +
    (geometry.ramaOutlierFraction == null ? 0 : Math.min(geometry.ramaOutlierFraction / 0.3, 1) * 0.15) +
    Math.min(geometry.chainBreaks * 0.08, 0.25);
  return clamp(1 - penalty, 0, 1);
}

function estimateAccuracy(
  geometry: Geometry,
  plddtMean: number | null,
  pae: PaeSummary,
): { point: number; low: number; high: number; domainGap: boolean } {
  let base = plddtMean == null ? 0.62 : clamp(plddtMean / 100, 0, 1);
  if (pae.longMean != null) base = 0.55 * base + 0.45 * clamp(1 - pae.longMean / 25, 0, 1);
  const clashRate = geometry.clashes / Math.max(geometry.residueCount, 1);
  base -= Math.min(clashRate / 0.2, 0.45);
  if (geometry.bondOutlierFraction != null) base -= Math.min(geometry.bondOutlierFraction, 0.25);
  base -= Math.min(geometry.chainBreaks * 0.03, 0.12);
  let point = clamp(base, 0.02, 0.98);
  let low = clamp(point - (plddtMean == null && pae.mean == null ? 0.16 : 0.08), 0, 1);
  let high = clamp(point + (plddtMean == null && pae.mean == null ? 0.16 : 0.08), 0, 1);
  const domainGap =
    pae.longMean != null && pae.localMean != null && pae.longMean - pae.localMean > 8 && pae.longMean > 15;
  if (domainGap) {
    point = Math.min(point, 0.6);
    low = Math.min(low, 0.42);
    high = Math.min(high, 0.74);
  }
  if (high < point) high = point;
  if (low > point) low = point;
  return { point, low, high, domainGap };
}

function pushGeometryReasons(reasons: string[], geometry: Geometry, residues: number): void {
  if (geometry.clashes === 0) reasons.push("No hard steric clashes.");
  else reasons.push(`${geometry.clashes} hard steric clash${geometry.clashes === 1 ? "" : "es"} across ${residues} residues.`);
  if (geometry.bondOutlierFraction != null && geometry.bondOutlierFraction > 0.02) {
    reasons.push(`${percent(geometry.bondOutlierFraction)} of covalent bonds are far from ideal length.`);
  }
  if (geometry.ramaOutlierFraction != null && geometry.ramaOutlierFraction > 0.1) {
    reasons.push(`${percent(geometry.ramaOutlierFraction)} of scored residues sit outside favored Ramachandran regions.`);
  }
  if (geometry.chainBreaks > 0) {
    reasons.push(
      `${geometry.chainBreaks} chain break${geometry.chainBreaks === 1 ? "" : "s"}, where neighboring Cα atoms are not about 3.8 Å apart.`,
    );
  }
}

function bondChecks(residues: Residue[]): boolean[] {
  const outliers: boolean[] = [];
  for (const residue of residues) {
    for (const [key, ideal] of Object.entries(BOND_IDEAL)) {
      const [left, right] = key.split("|");
      const distance = pairDistance(residue, residue, left, right);
      if (distance == null) continue;
      outliers.push(Math.abs(distance - ideal) > BOND_TOLERANCE);
    }
  }
  for (const [first, second] of consecutive(residues)) {
    const distance = pairDistance(first, second, "C", "N");
    if (distance == null) continue;
    outliers.push(Math.abs(distance - PEPTIDE_BOND) > BOND_TOLERANCE);
  }
  return outliers;
}

function angleChecks(residues: Residue[]): boolean[] {
  const outliers: boolean[] = [];
  for (const residue of residues) {
    for (const [key, ideal] of Object.entries(ANGLE_IDEAL)) {
      const [a, b, c] = key.split("|");
      const degrees = angleAt(residue, residue, residue, a, b, c);
      if (degrees == null) continue;
      outliers.push(Math.abs(degrees - ideal) > ANGLE_TOLERANCE);
    }
  }
  for (const [first, second] of consecutive(residues)) {
    const peptide = angleAt(first, first, second, "CA", "C", "N");
    if (peptide != null) outliers.push(Math.abs(peptide - PEPTIDE_ANGLE_CA_C_N) > ANGLE_TOLERANCE);
    const incoming = angleAt(first, second, second, "C", "N", "CA");
    if (incoming != null) outliers.push(Math.abs(incoming - PEPTIDE_ANGLE_C_N_CA) > ANGLE_TOLERANCE);
  }
  return outliers;
}

function ramaChecks(residues: Residue[]): boolean[] {
  const outliers: boolean[] = [];
  for (let index = 1; index < residues.length - 1; index += 1) {
    const previous = residues[index - 1];
    const residue = residues[index];
    const following = residues[index + 1];
    if (!bonded(previous, residue) || !bonded(residue, following)) continue;
    const prevC = previous.atoms.C;
    const nitrogen = residue.atoms.N;
    const ca = residue.atoms.CA;
    const carbon = residue.atoms.C;
    const nextN = following.atoms.N;
    if (!prevC || !nitrogen || !ca || !carbon || !nextN) continue;
    const phi = dihedral(prevC.xyz, nitrogen.xyz, ca.xyz, carbon.xyz);
    const psi = dihedral(nitrogen.xyz, ca.xyz, carbon.xyz, nextN.xyz);
    outliers.push(!ramaFavored(phi, psi, residue.resname));
  }
  return outliers;
}

function caTrace(residues: Residue[]): { breaks: number } {
  let breaks = 0;
  for (const [first, second] of consecutive(residues)) {
    const left = first.atoms.CA;
    const right = second.atoms.CA;
    if (!left || !right) continue;
    const distance = norm(sub(left.xyz, right.xyz));
    if (distance < CHAIN_BREAK_MIN || distance > CHAIN_BREAK_MAX) breaks += 1;
  }
  return { breaks };
}

function countClashes(residues: Residue[]): number {
  const coords: Vec3[] = [];
  const owners: number[] = [];
  let caOnly = 0;
  residues.forEach((residue, index) => {
    const heavy = Object.values(residue.atoms).filter((atom) => atom.element !== "H" && atom.element !== "D");
    if (heavy.length <= 1) caOnly += 1;
    for (const atom of heavy) {
      coords.push(atom.xyz);
      owners.push(index);
    }
  });
  if (coords.length < 2) return 0;
  const threshold = caOnly >= 0.8 * residues.length ? CA_ONLY_CLASH : HEAVY_CLASH;
  return gridClashes(coords, owners, threshold);
}

function gridClashes(coords: Vec3[], owners: number[], threshold: number): number {
  const buckets = new Map<string, number[]>();
  coords.forEach((xyz, index) => {
    const key = `${Math.floor(xyz[0] / threshold)},${Math.floor(xyz[1] / threshold)},${Math.floor(xyz[2] / threshold)}`;
    const bucket = buckets.get(key);
    if (bucket) bucket.push(index);
    else buckets.set(key, [index]);
  });
  let clashes = 0;
  const limit = threshold * threshold;
  for (const [key, members] of buckets) {
    const [x, y, z] = key.split(",").map(Number);
    for (let dx = -1; dx <= 1; dx += 1) {
      for (let dy = -1; dy <= 1; dy += 1) {
        for (let dz = -1; dz <= 1; dz += 1) {
          const neighbors = buckets.get(`${x + dx},${y + dy},${z + dz}`);
          if (!neighbors) continue;
          for (const left of members) {
            for (const right of neighbors) {
              if (right <= left) continue;
              if (Math.abs(owners[left] - owners[right]) < 2) continue;
              const delta = sub(coords[left], coords[right]);
              if (dot(delta, delta) < limit) clashes += 1;
            }
          }
        }
      }
    }
  }
  return clashes;
}

function resolvePlddt(residues: Residue[], override: number[] | null): number[] | null {
  if (override && override.length === residues.length) return override;
  let constant = 0;
  const values = residues.map((residue) => {
    const factors = Object.values(residue.atoms).map((atom) => atom.bfactor);
    if (Math.max(...factors) - Math.min(...factors) <= 0.15) constant += 1;
    return residue.atoms.CA?.bfactor ?? factors[0] ?? 0;
  });
  if (constant < 0.9 * residues.length) return null;
  if (values.some((value) => value < 0 || value > 100)) return null;
  return values;
}

type PaeSummary = { mean: number | null; localMean: number | null; longMean: number | null };

function emptyPae(): PaeSummary {
  return { mean: null, localMean: null, longMean: null };
}

function summarizePae(pae: number[][] | null): PaeSummary {
  if (!pae || pae.length < 2 || pae.length !== pae[0].length) return emptyPae();
  let meanSum = 0;
  let meanCount = 0;
  let localSum = 0;
  let localCount = 0;
  let longSum = 0;
  let longCount = 0;
  const size = pae.length;
  const stride = size > 1500 ? Math.ceil(size / 1500) : 1;
  for (let i = 0; i < size; i += stride) {
    for (let j = 0; j < size; j += stride) {
      if (i === j) continue;
      const value = pae[i][j];
      if (!Number.isFinite(value)) continue;
      const separation = Math.abs(i - j);
      meanSum += value;
      meanCount += 1;
      if (separation <= 8) {
        localSum += value;
        localCount += 1;
      }
      if (separation >= 12) {
        longSum += value;
        longCount += 1;
      }
    }
  }
  return {
    mean: meanCount ? meanSum / meanCount : null,
    localMean: localCount ? localSum / localCount : null,
    longMean: longCount ? longSum / longCount : null,
  };
}

function parseStructure(name: string, text: string): Atom[] {
  const trimmed = text.trim();
  if (!trimmed) throw new Error("That file is empty.");
  if (trimmed.length > 30_000_000) throw new Error("That file is larger than Foldcheck will read.");
  const lower = name.toLowerCase();
  const atoms =
    lower.endsWith(".cif") || lower.endsWith(".mmcif") || trimmed.startsWith("data_")
      ? parseCif(trimmed)
      : parsePdb(trimmed);
  if (atoms.length > MAX_ATOMS) {
    throw new Error(`This file has ${atoms.length} atoms. Foldcheck reads up to ${MAX_ATOMS}.`);
  }
  if (!atoms.length) throw new Error("No protein atoms found. Use a PDB or mmCIF coordinate file.");
  return atoms;
}

function parsePdb(text: string): Atom[] {
  const atoms: Atom[] = [];
  let seenModel = false;
  for (const raw of text.split(/\r?\n/)) {
    const record = raw.slice(0, 6).trim().toUpperCase();
    if (record === "MODEL") {
      if (seenModel && atoms.length) break;
      seenModel = true;
      continue;
    }
    if (record === "ENDMDL" && atoms.length) break;
    if (record !== "ATOM" && record !== "HETATM") continue;
    if (raw.length < 54) continue;
    const alt = raw.slice(16, 17).trim();
    if (alt !== "" && alt !== "A") continue;
    const atomName = raw.slice(12, 16).trim().toUpperCase();
    const resname = raw.slice(17, 20).trim().toUpperCase() || "UNK";
    if (record === "HETATM" || resname === "HOH" || resname === "WAT") continue;
    const element = raw.length >= 78 ? raw.slice(76, 78).trim().toUpperCase() : "";
    if (isHydrogen(atomName, element)) continue;
    const x = Number(raw.slice(30, 38));
    const y = Number(raw.slice(38, 46));
    const z = Number(raw.slice(46, 54));
    const resseq = Number(raw.slice(22, 26));
    if (![x, y, z, resseq].every(Number.isFinite)) continue;
    const bfactor = raw.length >= 66 && raw.slice(60, 66).trim() ? Number(raw.slice(60, 66)) : 0;
    atoms.push({
      name: atomName,
      element: element || elementFromName(atomName),
      resname,
      chain: raw.slice(21, 22).trim() || "A",
      resseq,
      icode: raw.slice(26, 27).trim(),
      xyz: [x, y, z],
      bfactor: Number.isFinite(bfactor) ? bfactor : 0,
    });
  }
  return atoms;
}

function parseCif(text: string): Atom[] {
  const lines = text.split(/\r?\n/);
  for (let index = 0; index < lines.length; index += 1) {
    if (lines[index].trim() !== "loop_") continue;
    index += 1;
    const columns: string[] = [];
    while (index < lines.length && lines[index].trim().startsWith("_")) {
      columns.push(lines[index].trim());
      index += 1;
    }
    if (!columns[0]?.startsWith("_atom_site.")) continue;
    const rows: string[][] = [];
    let current: string[] = [];
    while (index < lines.length) {
      const stripped = lines[index].trim();
      if (!stripped || stripped.startsWith("#")) {
        if (current.length) break;
        index += 1;
        continue;
      }
      if (stripped.startsWith("loop_") || stripped.startsWith("data_") || stripped.startsWith("_")) break;
      current.push(...cifTokens(stripped));
      index += 1;
      if (current.length >= columns.length) {
        rows.push(current.slice(0, columns.length));
        current = current.slice(columns.length);
      }
    }
    return atomsFromCif(columns, rows);
  }
  return [];
}

function atomsFromCif(columns: string[], rows: string[][]): Atom[] {
  const col = (...names: string[]) => {
    for (const name of names) {
      const token = name.startsWith("_") ? name : `_atom_site.${name}`;
      const found = columns.indexOf(token);
      if (found >= 0) return found;
    }
    return -1;
  };
  const groupI = col("group_PDB");
  const atomI = col("label_atom_id", "auth_atom_id");
  const elementI = col("type_symbol");
  const resnameI = col("label_comp_id", "auth_comp_id");
  const chainI = col("auth_asym_id", "label_asym_id");
  const seqI = col("auth_seq_id", "label_seq_id");
  const icodeI = col("pdbx_PDB_ins_code");
  const altI = col("label_alt_id");
  const xI = col("Cartn_x");
  const yI = col("Cartn_y");
  const zI = col("Cartn_z");
  const bI = col("B_iso_or_equiv");
  if (atomI < 0 || xI < 0 || yI < 0 || zI < 0) return [];
  const atoms: Atom[] = [];
  for (const row of rows) {
    if (groupI >= 0 && !["ATOM", "."].includes(row[groupI].toUpperCase())) continue;
    const name = row[atomI].trim().toUpperCase();
    let element = elementI >= 0 ? row[elementI].trim().toUpperCase() : "";
    if (element === "." || element === "?") element = "";
    if (isHydrogen(name, element)) continue;
    let resname = resnameI >= 0 ? row[resnameI].trim().toUpperCase() : "UNK";
    if (resname === "." || resname === "?") resname = "UNK";
    if (resname === "HOH" || resname === "WAT") continue;
    let alt = altI >= 0 ? row[altI].trim() : "";
    if (alt === "." || alt === "?") alt = "";
    if (alt !== "" && alt !== "A") continue;
    let chain = chainI >= 0 ? row[chainI].trim() : "A";
    if (!chain || chain === "." || chain === "?") chain = "A";
    let icode = icodeI >= 0 ? row[icodeI].trim() : "";
    if (icode === "." || icode === "?") icode = "";
    const seqToken = seqI >= 0 ? row[seqI].trim() : "1";
    if (seqToken === "." || seqToken === "?") continue;
    const x = Number(row[xI]);
    const y = Number(row[yI]);
    const z = Number(row[zI]);
    const resseq = Number(seqToken);
    if (![x, y, z, resseq].every(Number.isFinite)) continue;
    const bfactor = bI >= 0 && row[bI] !== "." && row[bI] !== "?" ? Number(row[bI]) : 0;
    atoms.push({
      name,
      element: element || elementFromName(name),
      resname,
      chain: chain.slice(0, 1) || "A",
      resseq,
      icode,
      xyz: [x, y, z],
      bfactor: Number.isFinite(bfactor) ? bfactor : 0,
    });
  }
  return atoms;
}

export function parseConfidence(text: string): { plddt: number[] | null; pae: number[][] | null } {
  let data: unknown;
  try {
    data = JSON.parse(text);
  } catch {
    throw new Error("The confidence file is not JSON.");
  }
  const plddt = findPlddt(data);
  const pae = findPae(data);
  if (!plddt && !pae) throw new Error("No pLDDT or PAE matrix was found in that JSON file.");
  return { plddt, pae };
}

function findPlddt(data: unknown): number[] | null {
  if (Array.isArray(data)) {
    if (data.length && isRecord(data[0])) return findPlddt(data[0]);
    return asVector(data);
  }
  if (!isRecord(data)) return null;
  for (const key of ["plddt", "pLDDT"]) {
    if (key in data) {
      const vector = asVector(data[key]);
      if (vector) return vector;
    }
  }
  for (const value of Object.values(data)) {
    const found = findPlddt(value);
    if (found) return found;
  }
  return null;
}

function findPae(data: unknown): number[][] | null {
  if (Array.isArray(data)) {
    if (data.length && isRecord(data[0])) {
      const record = data[0];
      if (Array.isArray(record.distance) && Array.isArray(record.residue1) && Array.isArray(record.residue2)) {
        const residue1 = record.residue1.map(Number);
        const residue2 = record.residue2.map(Number);
        const distance = record.distance.map(Number);
        const size = Math.max(...residue1, ...residue2);
        if (!Number.isFinite(size) || size < 2 || size > MAX_RESIDUES) return null;
        const matrix = Array.from({ length: size }, () => Array<number>(size).fill(Number.NaN));
        for (let index = 0; index < distance.length; index += 1) {
          const row = residue1[index] - 1;
          const column = residue2[index] - 1;
          if (row >= 0 && column >= 0 && row < size && column < size) matrix[row][column] = distance[index];
        }
        return matrix;
      }
      return findPae(record);
    }
    return asMatrix(data);
  }
  if (!isRecord(data)) return null;
  for (const key of ["predicted_aligned_error", "pae", "PAE"]) {
    if (key in data) {
      const matrix = asMatrix(data[key]);
      if (matrix) return matrix;
    }
  }
  return null;
}

function asVector(value: unknown): number[] | null {
  if (!Array.isArray(value) || value.length < 2) return null;
  const source = Array.isArray(value[0]) ? value[0] : value;
  if (!Array.isArray(source) || source.some((item) => typeof item === "object")) return null;
  const numbers = source.map(Number);
  if (numbers.length < 2 || numbers.some((item) => !Number.isFinite(item))) return null;
  return numbers;
}

function asMatrix(value: unknown): number[][] | null {
  if (!Array.isArray(value) || !value.length) return null;
  let rows = value;
  if (Array.isArray(rows[0]) && Array.isArray((rows[0] as unknown[])[0]) && Array.isArray((rows[0] as unknown[])[0])) {
    const nested = rows[0];
    if (Array.isArray(nested) && Array.isArray(nested[0])) rows = nested as unknown[];
  }
  if (!Array.isArray(rows[0]) || typeof (rows[0] as unknown[])[0] === "object") return null;
  const matrix = (rows as unknown[]).map((row) => (Array.isArray(row) ? row.map(Number) : []));
  const size = matrix.length;
  if (size < 2 || size > MAX_RESIDUES || matrix.some((row) => row.length !== size)) return null;
  if (matrix.some((row) => row.some((item) => !Number.isFinite(item)))) return null;
  return matrix;
}

export function buildPeptide(length: number, phi = HELIX[0], psi = HELIX[1]): Atom[] {
  let prevN: Vec3 = [0, 0, 0];
  let prevCa: Vec3 = [BOND_IDEAL["N|CA"], 0, 0];
  const dummy: Vec3 = [0, 1, 0];
  let prevC = placeAtom(dummy, prevN, prevCa, BOND_IDEAL["CA|C"], ANGLE_IDEAL["N|CA|C"], phi);
  const atoms: Atom[] = [];
  appendResidue(atoms, 1, prevN, prevCa, prevC);
  for (let resseq = 2; resseq <= length; resseq += 1) {
    const nextN = placeAtom(prevN, prevCa, prevC, PEPTIDE_BOND, PEPTIDE_ANGLE_CA_C_N, psi);
    const nextCa = placeAtom(prevCa, prevC, nextN, BOND_IDEAL["N|CA"], PEPTIDE_ANGLE_C_N_CA, 180);
    const nextC = placeAtom(prevC, nextN, nextCa, BOND_IDEAL["CA|C"], ANGLE_IDEAL["N|CA|C"], phi);
    prevN = nextN;
    prevCa = nextCa;
    prevC = nextC;
    appendResidue(atoms, resseq, prevN, prevCa, prevC);
  }
  return atoms;
}

function appendResidue(atoms: Atom[], resseq: number, nitrogen: Vec3, ca: Vec3, carbon: Vec3) {
  const oxygen = placeAtom(nitrogen, ca, carbon, BOND_IDEAL["C|O"], ANGLE_IDEAL["CA|C|O"], 180);
  const cb = placeAtom(carbon, nitrogen, ca, BOND_IDEAL["CA|CB"], CB_ANGLE, CB_DIHEDRAL);
  const records: [string, string, Vec3][] = [
    ["N", "N", nitrogen],
    ["CA", "C", ca],
    ["C", "C", carbon],
    ["O", "O", oxygen],
    ["CB", "C", cb],
  ];
  for (const [name, element, xyz] of records) {
    atoms.push({
      name,
      element,
      resname: "ALA",
      chain: "A",
      resseq,
      icode: "",
      xyz: [...xyz],
      bfactor: 0,
    });
  }
}

function assignBfactor(atoms: Atom[], value: number) {
  for (const atom of atoms) atom.bfactor = value;
}

function blockPae(length: number, local: number, cross: number): number[][] {
  const mid = length / 2;
  return Array.from({ length }, (_, row) =>
    Array.from({ length }, (_, column) => {
      if (row === column) return 0;
      const same = (row + 1 <= mid) === (column + 1 <= mid);
      return same ? local : cross;
    }),
  );
}

function placeAtom(a: Vec3, b: Vec3, c: Vec3, length: number, angleDeg: number, dihedralDeg: number): Vec3 {
  const angle = (angleDeg * Math.PI) / 180;
  const dihedral = (dihedralDeg * Math.PI) / 180;
  const bcHat = normalize(sub(c, b));
  const normal = cross(sub(b, a), sub(c, b));
  const nlen = norm(normal);
  if (nlen < 1e-8) throw new Error("Cannot place an atom from three colinear points");
  const nHat = scale(normal, 1 / nlen);
  const inPlane = cross(nHat, bcHat);
  return add(
    c,
    scale(bcHat, -length * Math.cos(angle)),
    scale(inPlane, length * Math.sin(angle) * Math.cos(dihedral)),
    scale(nHat, length * Math.sin(angle) * Math.sin(dihedral)),
  );
}

function groupResidues(atoms: Atom[]): Residue[] {
  const order: string[] = [];
  const groups = new Map<string, Residue>();
  for (const atom of atoms) {
    const key = `${atom.chain}|${atom.resseq}|${atom.icode}`;
    let residue = groups.get(key);
    if (!residue) {
      residue = { resname: atom.resname, chain: atom.chain, resseq: atom.resseq, icode: atom.icode, atoms: {} };
      groups.set(key, residue);
      order.push(key);
    }
    residue.atoms[atom.name] = atom;
  }
  return order.map((key) => groups.get(key)!);
}

function consecutive(residues: Residue[]): [Residue, Residue][] {
  const pairs: [Residue, Residue][] = [];
  for (let index = 1; index < residues.length; index += 1) {
    if (bonded(residues[index - 1], residues[index])) pairs.push([residues[index - 1], residues[index]]);
  }
  return pairs;
}

function bonded(first: Residue, second: Residue): boolean {
  return first.chain === second.chain && first.icode === second.icode && second.resseq - first.resseq === 1;
}

function pairDistance(left: Residue, right: Residue, leftName: string, rightName: string): number | null {
  const a = left.atoms[leftName];
  const b = right.atoms[rightName];
  if (!a || !b) return null;
  return norm(sub(a.xyz, b.xyz));
}

function angleAt(
  aRes: Residue,
  bRes: Residue,
  cRes: Residue,
  aName: string,
  bName: string,
  cName: string,
): number | null {
  const a = aRes.atoms[aName];
  const b = bRes.atoms[bName];
  const c = cRes.atoms[cName];
  if (!a || !b || !c) return null;
  const u = sub(a.xyz, b.xyz);
  const v = sub(c.xyz, b.xyz);
  const denom = norm(u) * norm(v);
  if (denom === 0) return null;
  return (Math.acos(clamp(dot(u, v) / denom, -1, 1)) * 180) / Math.PI;
}

function dihedral(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3): number {
  const b0 = sub(p0, p1);
  const b1raw = sub(p2, p1);
  const b2 = sub(p3, p2);
  const b1 = scale(b1raw, 1 / norm(b1raw));
  const v = sub(b0, scale(b1, dot(b0, b1)));
  const w = sub(b2, scale(b1, dot(b2, b1)));
  const x = dot(v, w);
  const y = dot(cross(b1, v), w);
  return (Math.atan2(y, x) * 180) / Math.PI;
}

function ramaFavored(phi: number, psi: number, resname: string): boolean {
  if (resname.toUpperCase() === "GLY") return true;
  if (phi >= -180 && phi <= -40 && ((psi >= 80 && psi <= 180) || (psi >= -180 && psi <= -160))) return true;
  if (phi >= -160 && phi <= -40 && psi >= -80 && psi <= 50) return true;
  if (phi >= 40 && phi <= 90 && psi >= -20 && psi <= 80) return true;
  return false;
}

function fraction(checks: boolean[]): number | null {
  if (!checks.length) return null;
  return checks.filter(Boolean).length / checks.length;
}

function elementFromName(name: string): string {
  const letters = name.toUpperCase().replace(/[^A-Z]/g, "");
  if (!letters) return "C";
  if (/^C[ABGDEZH]/.test(letters)) return "C";
  return letters[0];
}

function isHydrogen(name: string, element: string): boolean {
  return element === "H" || element === "D" || name.startsWith("H");
}

function cifTokens(line: string): string[] {
  const tokens: string[] = [];
  let index = 0;
  while (index < line.length) {
    if (/\s/.test(line[index])) {
      index += 1;
      continue;
    }
    if (line[index] === "'" || line[index] === '"') {
      const quote = line[index];
      index += 1;
      const start = index;
      while (index < line.length && line[index] !== quote) index += 1;
      tokens.push(line.slice(start, index));
      index += 1;
      continue;
    }
    const start = index;
    while (index < line.length && !/\s/.test(line[index])) index += 1;
    tokens.push(line.slice(start, index));
  }
  return tokens;
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === "object" && value !== null && !Array.isArray(value);
}

function average(values: number[]): number {
  return values.reduce((sum, value) => sum + value, 0) / values.length;
}

function percent(value: number): string {
  return `${Math.round(value * 100)}%`;
}

function clamp(value: number, low: number, high: number): number {
  return Math.min(high, Math.max(low, value));
}

function sub(a: Vec3, b: Vec3): Vec3 {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

function add(...vectors: Vec3[]): Vec3 {
  return vectors.reduce<Vec3>((sum, vector) => [sum[0] + vector[0], sum[1] + vector[1], sum[2] + vector[2]], [0, 0, 0]);
}

function scale(vector: Vec3, factor: number): Vec3 {
  return [vector[0] * factor, vector[1] * factor, vector[2] * factor];
}

function dot(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function cross(a: Vec3, b: Vec3): Vec3 {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

function norm(vector: Vec3): number {
  return Math.hypot(vector[0], vector[1], vector[2]);
}

function normalize(vector: Vec3): Vec3 {
  return scale(vector, 1 / norm(vector));
}
