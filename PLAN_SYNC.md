# PLAN_SYNC.md — Plan vs Reality Tracker

**Purpose:** Track any divergences between documentation (Memorial, SCRIPTS_PLAN, CHECKLIST,
RUNBOOK) and what is actually implemented. Reviewed and resolved at every sealed step.

**Rule:** Any divergence here must be communicated to the user before being resolved.
Resolving a divergence = update documentation OR update code, then delete the entry here.

---

## STATUS: CLEAN ✅

Last verified: 2026-02-28
Verified by: Claude (Option A cohort architecture approved by Rafael Botan)

No open divergences.

---

## Resolved divergences (archive)

### [RESOLVED 2026-02-28] GSE96058 treated as validation cohort — incorrect

- **Plan said:** SCAN-B (training) and GSE96058 (validation) split within the same GEO accession
- **Reality:** GSE96058 pData has NO training/validation split column. All 3069 samples are from SCAN-B. Decision: use ALL as training (Option A, approved by Rafael Botan 2026-02-28).
- **Files updated:** CLAUDE.md, MEMORY.md, Memorial_v6_1_CorePAM.md, COREPAM_REPRO_RUNBOOK_v2.md, CHECKLIST_MASTER_v6_1_CorePAM.md, SCRIPTS_PLAN_v6_1_CorePAM.md, PLAN_SYNC.md, METHODS_DRAFT.md, README.md
- **Deleted scripts:** 03_expression_preprocess_GSE96058.R, 06_zscore_and_score_GSE96058.R, 07_survival_analysis_GSE96058.R

---

### [RESOLVED 2026-02-28] Script naming: 01_ingest_raw_<COHORT>.R vs 01_download_raw_data.R

- **Plan said:** `01_ingest_raw_<COHORT>.R` (one script per cohort)
- **Reality:** `01_download_raw_data.R` (single monolithic script, all 4 OS cohorts)
- **Decision:** Update plan docs to match reality (Option A, approved by user)
- **Files updated:** `SCRIPTS_PLAN_v6_1_CorePAM.md`, `CHECKLIST_MASTER_v6_1_CorePAM.md`
- **Commit:** `0fb84a5`

---

## How to add a new divergence

When a discrepancy is detected, append a section here BEFORE fixing:

```markdown
### [OPEN yyyy-mm-dd] Brief description

- **Plan says:** ...
- **Reality:** ...
- **Impact:** ...
- **Awaiting:** user decision
```

After resolution, move to the "Resolved" archive above.

---

<!-- Entries are managed manually during review steps -->
