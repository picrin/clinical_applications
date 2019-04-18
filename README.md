This tool is implemented as a jupyter notebook and can be accessed at `http://dmbdi.adamkurkiewicz.com`. The site is password protected, to be granted access, please send an e-mail to `adam.kurkiewicz@glasgow.ac.uk`.

The template notebook is where you need to start, `DMBDI_explorer_template.ipynb`.

The notebook consists of cells. Each cell is numbered on the left (with numbering corresponding to the order of execution) and can be selected using a keyboard or a pointer. In order to execute a cell and select the next cell you can press shift-enter. The content of the cell consists of arbitrary python code, and can produce a text or a graphic output (e.g. an interactive plot).

The functionality of the tool I anticipate to be most useful is the ability to make interactive "railway" plots of gene expression, which should be useful for analysis of both APA and AS events.

That's how you can do it for muscle:

`plot_gene("TNNI1", muscle_metadata, "analysis/experiment_muscle", muscle_metadata.modal_allele)`

And that's how you can do it for blood:

`plot_gene("TNNI1", blood_metadata, "analysis/experiment_blood", blood_metadata.modal_allele)`

You can replace `"TNNI1"` with any valid HGNC gene name.

In principle the capabilities of the tool are limited only by the programming capabilities of its operator, and the underlying dataset. In practice, the following are supported:

1. Railway plots (analogous to Manhattan plots, but with microarray experiments).
2. Regressing expression data against:
   - Modal allele length from blood (e.g. `muscle_metadata.modal_allele`).
   - Statistically modelled progenitor allele length from blood (e.g. `muscle_metadata.progenitor_allele`).
   - Muscular Impairment Rating Scale (MIRS).
3. Load lists of genes previously implicated in DM1 for subsequent analysis with our DMBDI dataset.

We've also used the tool to build a pipeline which allows to evaluate efficacy of a potential DM1 treatment in a simulated clinical trial. This is described in README-analysis.md