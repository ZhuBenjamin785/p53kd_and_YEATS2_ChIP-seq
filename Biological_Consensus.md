# Biological consensus

## Bottom line

The current data support a real, perturbation-associated redistribution of H4K16ac, especially after p53 knockdown. They do **not yet show that p53 or YEATS2 changes 3D genome architecture**, and they do not support a shared p53–YEATS2 target set.

Hi-C is justified only if the next question is explicitly: **does either perturbation change chromatin contacts, compartments, TADs, or loops?** It is not yet justified as the automatic next step for explaining the H4K16ac peaks.

## Biological Consensus

In this system, p53 and YEATS2 appear to influence the H4K16ac chromatin state, but they do so through largely distinct genomic programs.

1. **p53 knockdown has the strongest chromatin effect.** The DiffBind summary reports 217 significant H4K16ac peaks at FDR < 0.05: 208 losses and 9 gains. The annotated files contain 208 losses and 5 gains, suggesting a small annotation/contig discrepancy that should be reconciled before final report. The direction is nevertheless unambiguous: loss dominates gain.

   <img width="436.8" height="316.8" alt="volcano_plot" src="https://github.com/user-attachments/assets/f03ce517-f5d4-4d5d-832c-695f36402810" />
   <img width="436.8" height="316.8" alt="MA_plot" src="https://github.com/user-attachments/assets/02c8888a-2320-4b9a-8316-502da1c8ddbe" />
   <img width="436.8" height="316.8" alt="peak_intensity_log2FC_distribution" src="https://github.com/user-attachments/assets/59f9ba61-39d8-4ef9-b9d8-6f7c39e05e5d" />



3. **The p53KD loss is largely promoter-proximal.** Of the 208 annotated loss peaks, 153 are within 1 kb of a TSS, with additional promoter peaks in the 1–2 kb and 2–3 kb windows. This is consistent with p53 knockdown changing H4K16ac at many promoter-associated regulatory regions. It does not by itself prove that transcription changed or that p53 directly binds every affected promoter.
   <img width="1312" height="792" alt="peak_annotation_summary_barplot" src="https://github.com/user-attachments/assets/0839b95a-e223-4e96-a6de-d5a93e79541c" />


5. **YEATS2 knockdown also changes H4K16ac, but on a smaller and different set of loci.** It has 61 significant peaks: 42 losses and 19 gains. Its significant peaks are distributed across intronic, promoter, distal-intergenic, exon, and UTR categories rather than showing the same strongly promoter-focused pattern.

6. **The p53KD and YEATS2KD effects are not the same local program.** The supplied overlap analysis finds 0 exact overlapping peak pairs among 217 p53KD and 61 YEATS2KD significant peaks. It also finds no shared nearby genes within 10 kb and no shared promoter-neighborhood genes. The one shared ChIPseeker-annotated gene is too little to establish biological convergence, and GO testing correctly reports the shared set as too small.

7. **The safest consensus is therefore “distinct chromatin responses,” not “one common p53–YEATS2 pathway.”** A possible model is that p53 and YEATS2 affect H4K16ac through different locus sets or regulatory routes. They may still converge at a broader chromatin or cellular state, but that is a hypothesis, not a result of the current overlap analysis.

8. **Expression-level conclusions are not currently available from the included exploratory RNA-seq tables.** The RNA metadata lists one sample for each of four conditions, and the exploratory fold-change tables have no p-values or adjusted p-values. They can suggest candidate directions, but they cannot provide replicated differential-expression evidence or establish an H4K16ac-to-expression chain.

9. **The p53 ChIP result is useful context but not a complete mechanism.** The 0-hour p53 ChIP has two replicates with reported FRiP values of 0.22 and 0.17, but the read depths are quite unequal (about 34.7 million versus 9.1 million reads). It supports examining p53 occupancy, but it does not establish that p53 directly controls the H4K16ac changes after knockdown.

## Why Hi-C could be the right next experiment

Hi-C measures DNA regions that are physically near one another in the nucleus. The current H4K16ac data measure a chromatin mark at genomic locations; they cannot distinguish among:

- a local promoter or enhancer change;
- a change in enhancer–promoter contact;
- a larger A/B compartment shift;
- a TAD-boundary or loop change; or
- a secondary consequence of altered transcription or cell state.

Hi-C would directly test the missing architectural layer. It becomes especially valuable if the intended model is that p53 or YEATS2 redistributes H4K16ac by changing enhancer–promoter communication or broader chromatin organization.

## Why moving directly to Hi-C may be premature

The current results do not contain a specific 3D-genome signature. Zero overlap between the p53KD and YEATS2KD peak sets argues against a simple shared local mechanism, but it does not predict a global architecture change. A Hi-C experiment is also expensive and can be difficult to interpret when there is only one sample per condition or when perturbations alter proliferation, cell composition, or other global states.

Hi-C would not tell us automatically which H4K16ac peak regulates which gene. Nor would a contact change prove that p53 or YEATS2 caused it. Without matched expression and chromatin-state evidence, a contact map could produce an interesting but weakly anchored result.

## Decision

**Recommendation: conditional go, but do not jump straight to a large Hi-C experiment yet.**

Move toward 3D-genome analysis if the project’s central question is architectural and the experiment can include:

- matched biological replicates for every perturbation and control;
- the same cell system, knockdown timing, and treatment state across assays;
- a design powered for the intended resolution, with biological QC before interpretation;
- RNA-seq with true replication, or another transcriptional readout;
- an orthogonal chromatin-state assay such as ATAC-seq or validated H4K16ac ChIP-seq; and
- targeted validation of a few candidate contacts, preferably at loci where H4K16ac, p53/YEATS2 occupancy, and expression point in the same direction.

If the immediate goal is simply to explain the current H4K16ac findings, the better next step is to validate the p53KD/YEATS2KD contrasts and connect them to replicated transcriptional changes. Then use targeted 3C/4C or a focused contact assay at selected loci before committing to genome-wide Hi-C.

## Claims to use and claims to avoid

Use:

> p53 and YEATS2 perturbation are associated with distinct H4K16ac redistribution patterns, with p53 knockdown showing a stronger promoter-proximal loss-dominant response.

Avoid:

> p53 and YEATS2 jointly remodel the 3D genome.

That stronger claim requires contact data and is not established by H4K16ac ChIP-seq, peak overlap, or the current exploratory RNA tables.

## Reconciliation items before finalizing the story

- Resolve why the p53KD summary reports 9 gains while the annotated gain table contains 5.
- Confirm the exact biological replicate structure and normalization used for every DiffBind contrast.
- Re-run or verify replicated RNA-seq statistics before making expression claims.
- Treat GO results based on two genes or very small gene sets as hypothesis-generating only.
- Keep spike-in-normalized ChIP signal distinct from matched-control log2(ChIP/Input) signal when describing quantitative effects.
