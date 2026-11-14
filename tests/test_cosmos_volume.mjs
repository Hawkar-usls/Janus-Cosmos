import test from "node:test";
import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {makeIllustration,parseFermiCatalog,fermiPoints,parseHI,bounds,scientificLimits} from "../docs/galaxy/model.mjs";

const original=JSON.parse(readFileSync(new URL("../data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-HE-PHOTONS-2026-10-18.json",import.meta.url),"utf8"));

test("gas model is deterministic, finite and explicitly illustrative",()=>{
  const a=makeIllustration({gasCount:150,starCount:20,seed:7});
  const b=makeIllustration({gasCount:150,starCount:20,seed:7});
  assert.deepEqual(a,b);
  assert.equal(bounds(a).count,170);
  assert.equal(bounds(a).finite,true);
  assert.equal(scientificLimits.galaxy,"ILLUSTRATIVE_DENSITY_MODEL_NOT_OBSERVATION");
  assert.ok(a.every(p=>p.kind==="model_gas"||p.kind==="model_star"));
});
test("Fermi source has 11 bursts and 38 actual events — none are sky positions",()=>{
  const events=parseFermiCatalog(original),pts=fermiPoints(events);
  assert.equal(original.bursts.length,11);assert.equal(events.length,38);
  assert.equal(new Set(events.map(e=>e.id)).size,38);
  assert.equal(Math.max(...events.map(e=>e.energy)),94.1);
  assert.equal(Math.min(...events.map(e=>e.energy)),10.4);
  assert.equal(Math.max(...events.map(e=>e.time)),34366);
  assert.equal(bounds(pts).finite,true);
  assert.ok(events.every(e=>e.observed&&!Object.hasOwn(e,"ra_deg")&&!Object.hasOwn(e,"distance_m")));
  assert.equal(scientificLimits.fermi,"OBSERVATIONS_IN_FEATURE_SPACE_NOT_ASTROMETRIC_POSITIONS");
});
test("schema and numerical faults do not silently produce scientific observations",()=>{
  assert.throws(()=>parseFermiCatalog({schema:"wrong",bursts:[]}));
  assert.throws(()=>parseFermiCatalog({schema:original.schema,bursts:[{grb:"TEST",redshift:"NaN",events:[]}]}));
  assert.throws(()=>makeIllustration({gasCount:90000}));
});
test("HI accepts ONLY explicit angular measurements, never inferred distance",()=>{
  const good={schema:"janus.cosmos.hi.angular_samples.v1",source_url:"https://example.org/hi4pi.fits",
    nhi_unit:"cm-2",samples:[{l_deg:5,b_deg:3,nhi_cm2:1e21}]};
  const hi=parseHI(good);assert.equal(hi.length,1);
  assert.ok(Math.abs(Math.hypot(hi[0].x,hi[0].y,hi[0].z)-4.25)<1e-9);
  assert.equal(hi[0].kind,"hi_angular");
  assert.ok(!Object.hasOwn(hi[0],"distance_pc"));
  assert.throws(()=>parseHI({...good,nhi_unit:"K"}));
  assert.throws(()=>parseHI({...good,samples:[{l_deg:361,b_deg:2,nhi_cm2:1e21}]}));
  assert.throws(()=>parseHI({...good,source_url:"javascript:alert(1)"}));
});
test("viewer remains additive and does not touch original observer",()=>{
  const html=readFileSync(new URL("../docs/galaxy/index.html",import.meta.url),"utf8");
  const view=readFileSync(new URL("../docs/galaxy/view.mjs",import.meta.url),"utf8");
  assert.match(html,/JANUS COSMOS/);assert.match(html,/canvas id="cosmos"/);
  assert.match(view,/getContext\("webgl2"/);
  assert.match(view,/parseFermiCatalog/);assert.match(view,/parseHI/);
  assert.match(view,/FERMI.*JSON|Fermi.*JSON/);
  assert.doesNotMatch(view,/Math\.random\(/);
});
