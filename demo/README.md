# IBMT Clinical — demos

Two working prototypes for running clinical trials in Africa, plus the model code behind
the first one.

| | What it does | Path |
|---|---|---|
| **Pharmacogenomic Dose Simulator** | Predicts how a standard dose lands in West African vs European patients, from published pharmacogenomic data | `simulator/` |
| **Regulatory Drafting** | Turns a sponsor protocol into a Nigerian submission package, every requirement traced to a public source | `regulatory/` |
| **Model** | The fit and validation behind the simulator | `model/model.py` |

Both demos are single self-contained HTML files with no build step and no dependencies.

---

## Running the model

```bash
pip install -r model/requirements.txt
python3 model/model.py
```

It prints every number shown on the simulator page, including the fitted constants that the
page hard-codes. If the page and the model ever disagree, the model is the source of truth.

## What the simulator claims, and what actually supports it

Two gene–drug pairs, chosen because the effect runs in **opposite directions**:

- **Efavirenz / CYP2B6** — West African patients are *over*-exposed at the standard 600 mg,
  because slow-metabolising genotypes are far more common. About one Nigerian in five is a
  slow or ultra-slow metaboliser, against roughly one in sixteen of European ancestry.
- **Tacrolimus / CYP3A5** — the reverse. Expressers clear the drug faster and need *more*,
  and nearly every West African patient is an expresser while few European patients are.

For each drug the model is fitted on published concentration data only. Dosing guidelines
(CPIC 2019 for efavirenz, CPIC 2015 for tacrolimus) were **not used in fitting** and serve as
held-out validation. The calibration table on the page bands the two groups separately,
because rows used to set the parameters cannot test the model — they match by construction.

The strongest single result is the efavirenz GT/GG clearance ratio: fitted on Nigerian adult
data, it predicts 0.714 against 0.706 observed in a different cohort of Zambian and Ugandan
*children*. Different region, different age group, 1.1% apart.

## Known limitations

Stated here rather than buried, because they are the honest boundary of what a public-data
model can do.

**Efavirenz**
- Only two CYP2B6 loci are modelled; other variants affect clearance.
- The 983T>C clearance values are scaled from a paediatric population PK study, not fitted
  on adults.
- Indexed on the 983 locus rather than resolving the full CYP2B6\*18 haplotype. In the same
  West African sample the \*18 haplotype frequency was 4.2% against 4.7% for the 983C allele,
  so the two are close but not identical.

**Tacrolimus**
- The 10–15 ng/mL target is narrow against ~70% coefficient of variation, so in-target
  fractions are low for every cohort at every dose. This is why tacrolimus requires
  therapeutic drug monitoring, and why the useful output is the *difference between cohorts*
  rather than absolute target attainment.
- The lognormal approximation runs about 6–9% high at the upper quartile, because the
  published distribution's median sits above its geometric midpoint.
- CPIC notes that people of sub-Saharan African descent need high tacrolimus doses
  regardless of CYP3A5 genotype, so factors beyond this gene are at work and a genotype-only
  model under-predicts African requirements. The data needed to explain the rest has never
  been generated — which is the gap this company exists to close.

**Populations.** There is no "African" population in the simulator, by design. How much
allele frequency varies across Africa depends on the gene: for CYP3A5 the spread is
twentyfold, and Ethiopia sits closer to European ancestry than to Nigeria. Country labels
name the country a sampled cohort came from, not the whole country. Where a frequency has
not been published — Ghana and Senegal for CYP2B6 516G>T — the population is listed but not
offered for simulation, rather than filled in by assumption.

---

**Not regulatory or medical advice.** These are prototypes built to test an approach.
