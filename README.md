# vitamin-a

This is code and paper repository for the work-in-progress paper on vitamin A supplementation and child mortality.

The paper is in the main folder

Code structure is very simple, there are three parts:

1. `simulation` (outputs in `simulation_results/`): simulations comparing performance of FE and RE models; this is not yet fully clean yet, but the code should run stand-alone in R
2. `meta-analysis`: the main script that replicates all analysis is `fit_ma.R`
3. `paper`: this has scripts used to generate tables, figures and some text (not fully cleaned yet)

## Data

In `input/` folder I included two files:
- `imdad2022.csv` is a manual transcription of risk ratios in Figure 3 in Imdad et al 2022
- `mayo_wilson_fig3.xlsx` is the same, done for rate ratios in Figure 3 in Mayo-Wilson et al 2011



## References

Imdad, Aamer, Evan Mayo-Wilson, Maya R. Haykal, Allison Regan, Jasleen Sidhu, Abigail Smith, and Zulfiqar A. Bhutta. ‘Vitamin A Supplementation for Preventing Morbidity and Mortality in Children from Six Months to Five Years of Age’. Cochrane Database of Systematic Reviews, no. 3 (2022). https://doi.org/10.1002/14651858.CD008524.pub4.

Mayo-Wilson, Evan, Aamer Imdad, Kurt Herzer, Mohammad Yawar Yakoob, and Zulfiqar A. Bhutta. ‘Vitamin A Supplements for Preventing Mortality, Illness, and Blindness in Children Aged under 5: Systematic Review and Meta-Analysis’. BMJ (Clinical Research Ed.) 343 (25 August 2011): d5094. https://doi.org/10.1136/bmj.d5094.
