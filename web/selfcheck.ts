import { scoreSample, scoreText } from "./score.ts";

const helix = scoreSample("helix");
const clash = scoreSample("clash");
const shift = scoreSample("shift");

function assert(condition: boolean, message: string) {
  if (!condition) throw new Error(message);
}

assert(helix.verdict === "reliable", `helix verdict ${helix.verdict} viability ${helix.viability}`);
assert(helix.clashes === 0, `helix clashes ${helix.clashes}`);
assert(helix.chainBreaks === 0, `helix breaks ${helix.chainBreaks}`);
assert((helix.bondOutlierFraction ?? 0) === 0, `helix bonds ${helix.bondOutlierFraction}`);
assert(clash.verdict === "reject", `clash verdict ${clash.verdict}`);
assert(clash.clashes > 0, "clash sample produced no clashes");
assert((clash.plddtMean ?? 0) > 80, "clash sample lost its high pLDDT");
assert(clash.viability < helix.viability, "clash scored as viable as the helix");
assert(shift.verdict !== "reliable", `shift looked reliable ${shift.predictedLddt}`);
assert(shift.predictedLddt < helix.predictedLddt - 0.15, `shift lDDT ${shift.predictedLddt}`);
assert(shift.chainBreaks >= 1, "shift did not break the chain");

const pdb = `ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       2.000   1.200   0.000  1.00 90.00           C
ATOM      4  N   ALA A   2       3.200   1.400   0.500  1.00 40.00           N
ATOM      5  CA  ALA A   2       4.500   1.600   0.800  1.00 40.00           C
ATOM      6  C   ALA A   2       5.200   2.800   0.700  1.00 40.00           C
END
`;
const parsed = scoreText("tiny.pdb", pdb, null);
assert(parsed.residues === 2, `parsed ${parsed.residues} residues`);
assert(parsed.plddtMean != null && parsed.plddtMean < 70, "pLDDT was not read from B-factors");

const pae = JSON.stringify({ predicted_aligned_error: [[0, 4], [4, 0]], plddt: [91, 88] });
const withPae = scoreText("tiny.pdb", pdb, pae);
assert(withPae.paeMean != null, "PAE JSON was ignored");
assert((withPae.plddtMean ?? 0) > 80, "JSON pLDDT did not override B-factors");

console.log(
  JSON.stringify(
    {
      helix: [helix.verdict, helix.viability, helix.predictedLddt, helix.geometryScore],
      clash: [clash.verdict, clash.viability, clash.clashes, clash.plddtMean],
      shift: [shift.verdict, shift.predictedLddt, shift.paeLong, shift.chainBreaks],
    },
    null,
    2,
  ),
);
