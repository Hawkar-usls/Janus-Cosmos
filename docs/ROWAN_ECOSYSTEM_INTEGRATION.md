# JANUS Cosmos ↔ Rowan / TOPA / AURA / Demiurge

**Role:** cross-domain *evidence transport and reproducibility*, not transfer
of biomedical findings into astronomical or geophysical evidence.

The independent Cosmos witness reads seven public JSON documents from three
repositories plus the Meta Registry: Rowan canonical progress, AURA hypotheses,
Demiurge council, coordination latest state and its immutable receipt, TOPA
sequence-shadow independent replay, and the pinned 44-feature shadow model.

It validates the exact canonical-JSON SHA-256 links used by the source
protocols. TOPA's model-file checksum is instead SHA-256 of the original
model file bytes; these hash schemes must not be confused.

The source audit checks:

- Rowan run, 8-case progress and upstream artifact bindings all agree.
- AURA and Demiurge packets match Meta Registry by content hash.
- Coordination latest pointer matches its immutable source-binding receipt.
- The independent TOPA replay corresponds to the archived sequence-shadow
  model, and remains **SHADOW_ONLY** with external validation open.
- No component authorizes new studies, Rowan provider jobs, paid credits,
  production model promotion, or empirical proof from proposed hypotheses.

A passing check writes a deterministic receipt to
`data/integrations/rowan/receipts/COSMOS-ROWAN-<run>-<key>.json`
and `data/integrations/rowan/LATEST.json`. Existing immutable receipts cannot
be silently rewritten. Replaying the same inputs does not change the repo.
If upstream binding or source access fails, the workflow fails **closed**.

GitHub Actions `cosmos-rowan-ecosystem-witness.yml` runs every six hours,
on manual dispatch, and on bridge code changes. On PRs it is **read-only**.
On default-branch runs only a separate tightly-scoped writer may archive
receipts. There is no Rowan API key, provider submission, or expenditure.

For local checks:

```bash
python -m unittest discover -v -s tests -p 'test_rowan_cosmos_bridge.py'
python tools/rowan_cosmos_bridge.py --out /tmp/cosmos-rowan-audit.json
```

**Isolation:** this integration does not modify `PROJECT_STATUS.json`,
`CANONICAL_RECEIPT_POINTER.json`, ocean/Cousteau branches, Earth/Space bridge,
OSIRIS proof verdicts, BHB/COSMOS astrophysics gates, or the scientific
authority of TOPA, AURA, or Demiurge. It may share source-verification and
negative-control *methods*; it cannot transform protein scores into evidence
for astronomical, extraordinary, or complexity-theory claims.
