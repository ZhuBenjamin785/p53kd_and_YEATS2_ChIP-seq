# Biological consensus

## Bottom line

The current data support a real, perturbation-associated redistribution of H4K16ac, especially after p53 knockdown. They do not yet show that p53 or YEATS2 changes 3D genome architecture, and they do not support a shared p53–YEATS2 target set.


## Biological Consensus

In this system, p53 and YEATS2 appear to influence the H4K16ac chromatin state, but they do so through largely distinct genomic programs.

1. P53 knockdown has the strongest effect on H4K16ac. The DiffBind Summary found 217 significant H4K16ac peaks above 0.05 FDR. the summary found 208 losses and 9 gains, however the annotated files found 208 losses and only 5 gains, which does suggest a discrepancy that should be reviewed, however the it is nevertheless clear that loss dominates gain in p53KD.
   <img width="436.8" height="316.8" alt="volcano_plot" src="https://github.com/user-attachments/assets/f03ce517-f5d4-4d5d-832c-695f36402810" />
   <img width="436.8" height="316.8" alt="MA_plot" src="https://github.com/user-attachments/assets/02c8888a-2320-4b9a-8316-502da1c8ddbe" />
   <img width="436.8" height="316.8" alt="peak_intensity_log2FC_distribution" src="https://github.com/user-attachments/assets/59f9ba61-39d8-4ef9-b9d8-6f7c39e05e5d" />



2. The discovered loss in p53KD is largely near promoter sites. Of the 208 loss peaks I found, 153 were found within 1kb of a Trascription Start Site (TSS) with additional loss peaks in that 1-2kb and 2-3kb window. This is consistent with the idea of p53 knockdown changing H4k16ac at many promoter related regulatory regions. It does not, however, prove that transcription changed or that p53 directly binds every affected promoter.
   <img width="437" height="264" alt="peak_annotation_summary_barplot" src="https://github.com/user-attachments/assets/0839b95a-e223-4e96-a6de-d5a93e79541c" />
   <img width="389" height="264" alt="TSS_average_profile" src="https://github.com/user-attachments/assets/f283ad92-c63a-45c1-9d0d-d21f263c1574" />




3. Knocking down YEATS2 was also found to change H4K16ac, however it was on smaller and different loci. I found 61 significant peaks, 42 were losses and 19 were gains, distributed more evenly across intronic, promoter, distal intergenic, exon and UTR, compared to p53 which showed the same strongly promoter focused pattern.
<img width="364" height="264" alt="peak_intensity_log2FC_distribution" src="https://github.com/user-attachments/assets/64f3b3a0-5ebd-4311-9a93-233a3f57e3b4" />
<img width="437" height="264" alt="peak_annotation_summary_barplot" src="https://github.com/user-attachments/assets/46ea9ab1-8392-4127-b288-9f84d5415c08" />
<img width="389" height="264" alt="TSS_average_profile" src="https://github.com/user-attachments/assets/22d3e698-62ae-4256-8749-0f056d4d2add" />
<img width="364" height="264" alt="volcano_plot" src="https://github.com/user-attachments/assets/6a1104af-6e51-4c88-bdf0-e0047617b801" />
<img width="364" height="264" alt="MA_plot" src="https://github.com/user-attachments/assets/9a506a07-4b04-458b-8049-4820cc7ee04d" />



4. I found the effects of p53KD and YEATS2KD to not be on the same local program. When an overlap analysis was conducted, 0 exact overlapping peak pairs were found amoung 217 p53KD and 61 YEATS2KD significant peaks. Additionally, it finds no shared nearby genes within 10 kb and no shared promoter-neighborhood genes. I did find 1 shared ChIPseeker-annotated gene, however it was too small to establish any type of biological convergence. GO testing also reported the shared set as too small.
<img width="352" height="297" alt="shared_gene_overlap_venn" src="https://github.com/user-attachments/assets/68363f55-ef6b-4d8a-9f29-d2274b10219b" />
<img width="352" height="297" alt="promoter_gene_overlap_venn" src="https://github.com/user-attachments/assets/30292f9c-a760-43de-adab-30eec1da9e79" />
<img width="352" height="297" alt="exact_peak_overlap_venn" src="https://github.com/user-attachments/assets/f21c5bfb-4d92-465a-93e1-7c800f464ee6" />


7. The safest consensus is therefore “distinct chromatin responses,” not “one common p53–YEATS2 pathway.” A possible model is that p53 and YEATS2 affect H4K16ac through different locus sets or regulatory routes. They may still converge at a broader chromatin or cellular state, but that is a hypothesis, not a result of the current overlap analysis.

8. Expression-level conclusions are not currently available from the included exploratory RNA-seq tables. The RNA metadata lists one sample for each of four conditions, and the exploratory fold-change tables have no p-values or adjusted p-values. They can suggest candidate directions, but they cannot provide replicated differential-expression evidence or establish an H4K16ac-to-expression chain.

9. The p53 ChIP result is useful context but not a complete mechanism. The 0-hour p53 ChIP has two replicates with reported FRiP values of 0.22 and 0.17, but the read depths are quite unequal (about 34.7 million versus 9.1 million reads). It supports examining p53 occupancy, but it does not establish that p53 directly controls the H4K16ac changes after knockdown.



