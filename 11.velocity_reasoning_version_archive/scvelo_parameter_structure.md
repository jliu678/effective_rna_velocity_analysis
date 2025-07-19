---
title: "Scvelo Parameter Structure"
date: 2025-07-10
draft: True
---

# A) scVelo parameter structure

## Gene-specific parameters (same for all cells)

* `.var['fit_alpha']` - Transcription rate for each gene ({{< math >}} $\alpha_g$ {{< /math >}})
* `.var['fit_beta']` - Splicing rate for each gene ({{< math >}} $\beta_g$ {{< /math >}})  
* `.var['fit_gamma']` - Degradation rate for each gene ({{< math >}} $\gamma_g$ {{< /math >}})
* `.var['fit_t_']` - Switching time for each gene ({{< math >}} $t_{switch,g}$ {{< /math >}})

**Dimensions**: {{< math >}} $(n_{genes},)$ {{< /math >}} - one value per gene

## Cell-and-gene-specific data

* `.layers['fit_t']` - Latent time varies by both cell AND gene ({{< math >}} $t_{c,g}$ {{< /math >}})
* `.layers['velocity']` - Velocity varies by both cell AND gene ({{< math >}} $v_{c,g}$ {{< /math >}})

**Dimensions**: {{< math >}} $(n_{cells} \times n_{genes})$ {{< /math >}} - different value for each cell-gene combination

## Mathematical representation

For a given gene {{< math >}} $g$ {{< /math >}} and cell {{< math >}} $c$ {{< /math >}}:

{{< math >}}
$$
\begin{align}
\alpha_g &= \text{transcription rate (gene-specific)} \\
\beta_g &= \text{splicing rate (gene-specific)} \\
\gamma_g &= \text{degradation rate (gene-specific)} \\
t_{c,g} &= \text{latent time (cell- and gene-specific)} \\
v_{c,g} &= \text{velocity (cell- and gene-specific)}
\end{align}
$$
{{< /math >}}

## Why this makes sense biologically

* Each gene has **intrinsic kinetic properties** - the rates {{< math >}} $\alpha_g$, $\beta_g$, $\gamma_g$ {{< /math >}} describe how fast the gene transcribes, splices, and degrades
* But different cells are at **different stages** of the process for each gene - the latent time {{< math >}} $t_{c,g}$ {{< /math >}} varies
* The **rates are gene properties**, the **timing/stage varies by cell**

This means:
- If you want to compare kinetic rates between genes → look at `.var`
- If you want to see where different cells are in the process for a specific gene → look at `.layers['fit_t']`

## Data access patterns

```python
# Same gene has same intrinsic kinetic rates across all cells
gene_alpha = adata.var['fit_alpha']['Ins1']  # α_Ins1 (same for all cells)
gene_beta = adata.var['fit_beta']['Ins1']    # β_Ins1 (same for all cells)  
gene_gamma = adata.var['fit_gamma']['Ins1']  # γ_Ins1 (same for all cells)

# But latent time varies per cell for the same gene
ins1_latent_times = adata.layers['fit_t'][:, gene_idx]  # t_c,Ins1 (different per cell)
# Cell 1 might be at t=0.2, Cell 2 at t=0.8 for the same gene
```


# B) Likelihood in scVelo

## Model Fit

The likelihood, specifically {{< math >}} $\text{fit\_likelihood}$ {{< /math >}} in scVelo, quantifies how well the observed unspliced and spliced mRNA counts for a given gene across all cells fit the gene's inferred kinetic model. This kinetic model describes the interplay of transcription, splicing, and degradation rates, and how these rates determine the accumulation and decay of unspliced and spliced mRNA over time.

## Expectation-Maximization (EM)

scVelo's dynamical model uses an EM framework. In this iterative process, the parameters of the splicing kinetics are estimated for each gene:
- Rate of transcription: {{< math >}} $\alpha$ {{< /math >}}
- Rate of splicing: {{< math >}} $\beta$ {{< /math >}}
- Rate of degradation: {{< math >}} $\gamma$ {{< /math >}}

Simultaneously, a cell-internal "latent time" (representing the progression along a developmental or biological process) is inferred for each cell. The likelihood is maximized during this process through the optimization:

{{< math >}} $$\max_{\alpha, \beta, \gamma, t} \mathcal{L}(\alpha, \beta, \gamma, t | u, s)$$ {{< /math >}}

where {{< math >}} $u$ {{< /math >}} represents unspliced counts and {{< math >}} $s$ {{< /math >}} represents spliced counts.

## Phase Portrait

Conceptually, scVelo tries to fit a "phase portrait" for each gene, which plots unspliced mRNA levels against spliced mRNA levels. The dynamic model aims to find a trajectory within this phase portrait that best explains the data points (cells). A high likelihood means the observed cells' mRNA levels for that gene fall very close to the trajectory predicted by the model.

The kinetic equations governing this relationship are:

{{< math >}} $$\frac{du}{dt} = \alpha(t) - \beta u$$ {{< /math >}}

{{< math >}} $$\frac{ds}{dt} = \beta u - \gamma s$$ {{< /math >}}

## Why do Genes with High Likelihood Matter for Trajectory Inference?

### Reliable Dynamics

Genes with high likelihood scores are those for which scVelo can confidently infer the underlying kinetic parameters and their dynamic behavior. If a gene has a low likelihood, it means the model struggles to explain its splicing dynamics, possibly due to:
- Noise in the data
- Violations of model assumptions
- Lack of clear dynamic changes

High-likelihood genes, therefore, offer more reliable insights into the actual cellular processes.

### Driving Force of Change (Driver Genes)

Genes that exhibit strong, consistent dynamic behavior (e.g., clear induction or repression phases) are often the "driver genes" of a developmental trajectory. These are the genes whose expression changes are actively pushing cells from one state to another. scVelo systematically detects these driver genes by their high likelihood in the dynamic model. Their clear kinetic profiles provide strong directional information.

### Accurate Latent Time and Velocity Estimation

The overall accuracy of the inferred developmental trajectory and RNA velocities relies on information from many genes. Genes with high likelihood contribute more robust and informative signals to the collective inference of:
- Cell-specific latent time {{< math >}} $t_c$ {{< /math >}}
- Directional RNA velocity vectors {{< math >}} $\vec{v}$ {{< /math >}}

If a gene's dynamics are poorly modeled (low likelihood), its contribution to the overall velocity field would be noisy or misleading.

### Biological Relevance

Biologically, high-likelihood genes are often those that are actively involved in the differentiation process being studied. They are likely to be:
- Transcription factors
- Signaling molecules
- Structural genes whose precise regulation is critical for cell fate transitions

By focusing on these genes, researchers can identify the key molecular players orchestrating cellular development.

## Practical Implications of the Threshold (>0.1)

The threshold of "{{< math >}} $> 0.1$ {{< /math >}}" (or any specific value) is typically an empirical cutoff. It means that genes with a {{< math >}} $\text{fit\_likelihood}$ {{< /math >}} score above this value are considered to have a sufficiently strong and well-explained dynamic signal to be used for further analysis.

### Genes Below This Threshold Might Be:

#### Non-dynamic
They might be expressed constitutively and not show significant changes over the developmental process. These genes have velocities {{< math >}} $v \approx 0$ {{< /math >}}.

#### Noisy
Their counts might be too low or too variable to accurately infer kinetic parameters. The signal-to-noise ratio is insufficient for reliable parameter estimation.

#### Poorly Modeled
Their kinetics might not fit the simple splicing model assumed by scVelo. The actual biological process may involve more complex regulatory mechanisms not captured by the basic kinetic equations.

## Result

By filtering for high-likelihood genes using the criterion:

{{< math >}} $$\text{selected genes} = \{g : \text{fit\_likelihood}_g > 0.1\}$$ {{< /math >}}

scVelo focuses its trajectory inference on the most informative and reliable genes, leading to more accurate and biologically meaningful results for understanding cell fate decisions and developmental paths.

This filtering approach ensures that the downstream analysis of cellular trajectories is based on genes with:
- Strong dynamic signals
- Reliable kinetic parameter estimates  
- Clear directional information for velocity computation
- Biological relevance to the differentiation process


# C) Apply Kinetic Parameters depending on Cell Type

Traditional RNA velocity (and even the basic scVelo dynamical model) often assumes that a gene's kinetic parameters ({{< math >}} $\alpha, \beta, \gamma$ {{< /math >}}) are constant across all cells, or at least across a continuous trajectory. However, this assumption can be violated in complex biological systems for several reasons:

## Distinct Cell Types/Lineages

Different cell types are fundamentally different. A fibroblast and a neuron, even if derived from a common progenitor, will have vastly different gene regulatory networks. Genes expressed in both might be regulated by different transcription factors, leading to different speeds of induction, splicing, or degradation.

## Governed by Different Network Structures

The set of transcription factors, chromatin modifiers, and RNA-binding proteins that control a gene's expression can change dramatically as a cell differentiates or switches lineage. These changes in the regulatory network directly impact the kinetic rates.

## Alternative Splicing and Polyadenylation

### Alternative Splicing
A single gene can produce multiple mRNA isoforms through alternative splicing. Different isoforms might have different splicing rates, degradation rates, or even different translation efficiencies. If different cell types preferentially use different isoforms, the "effective" kinetic rates of that gene would appear to change.

### Alternative Polyadenylation (APA)
Variations in polyadenylation sites can lead to mRNA transcripts with different 3' UTR lengths. Longer or shorter 3' UTRs can impact mRNA stability and degradation rates, leading to differential kinetics.

## Modulations in Degradation

mRNA degradation rates are not static. They can be actively regulated by microRNAs, RNA-binding proteins, and other cellular factors. If these regulatory elements are differentially expressed or active in different cell types, the degradation kinetics of target genes will vary.

## Compromised Model Fit

If a single set of kinetic parameters is forced onto a gene whose dynamics actually differ across cell types, the model's overall fit (likelihood) for that gene will be poor. This could obscure true dynamic information.

# The Solution: Differential Kinetic Test (Likelihood Ratio Test)

scVelo addresses this by performing a likelihood ratio test (LRT) for differential kinetics. This is a standard statistical approach for comparing nested models:

## The Null Hypothesis (Simpler Model)

For a given gene and a set of cell clusters/lineages, the null hypothesis is that all cells follow a single kinetic model for that gene. This means a single set of parameters {{< math >}} $(\alpha,\beta,\gamma)$ {{< /math >}} is sufficient to explain the unspliced and spliced mRNA levels across all cells in all groups. The maximum likelihood under this assumption is {{< math >}} $L(\hat{\theta}_0)$ {{< /math >}}, where {{< math >}} $\hat{\theta}_0$ {{< /math >}} represents the parameters fit to the entire population.

## The Alternative Hypothesis (More Complex Model)

The alternative hypothesis is that each distinct cluster/lineage exhibits its own independent kinetic model for that gene. This means a separate set of parameters {{< math >}} $(\alpha_k, \beta_k, \gamma_k)$ {{< /math >}} is fit for each cluster {{< math >}} $k$ {{< /math >}}. The maximum likelihood under this assumption is {{< math >}} $L(\hat{\theta})$ {{< /math >}}, where {{< math >}} $\hat{\theta}$ {{< /math >}} represents the parameters fit to each subgroup independently. This is a "more complex" model because it has more parameters (e.g., {{< math >}} $3 \times \text{number of clusters}$ {{< /math >}} parameters instead of just 3).

## The Likelihood Ratio (LR) Test Statistic

{{< math >}} $$\text{LR} = -2 \ln \frac{\sup_{\theta_0} L(\theta_0)}{\sup_{\theta} L(\theta)}$$ {{< /math >}}

Let's break this down:

- {{< math >}} $\sup_{\theta_0} L(\theta_0)$ {{< /math >}}: This is the maximum likelihood of the simpler model (null hypothesis). It's the highest probability of observing the data if a single kinetic model applies to all groups.

- {{< math >}} $\sup_{\theta} L(\theta)$ {{< /math >}}: This is the maximum likelihood of the more complex model (alternative hypothesis). It's the highest probability of observing the data if separate kinetic models apply to each group.

- **Ratio**: The ratio {{< math >}} $\frac{\sup_{\theta_0} L(\theta_0)}{\sup_{\theta} L(\theta)}$ {{< /math >}} will always be {{< math >}} $\leq 1$ {{< /math >}}. If the single model is a good fit, this ratio will be close to 1. If the independent fits are much better, this ratio will be small.

- **Log Transformation**: Taking the natural logarithm transforms the ratio into a difference.

- **Multiplication by -2**: The {{< math >}} $-2\ln(\text{ratio})$ {{< /math >}} ensures that the LR test statistic is positive and follows a chi-squared ({{< math >}} $\chi^2$ {{< /math >}}) distribution under the null hypothesis (with degrees of freedom equal to the difference in the number of parameters between the two models).

## Interpretation

### Large LR Value
A large LR value (and thus a small p-value) indicates that the more complex model (independent kinetics per cluster) provides a significantly better fit to the data than the simpler model (single kinetics for all clusters). This means that the gene indeed exhibits differential kinetics across the tested clusters.

### Small LR Value
A small LR value (large p-value) suggests that there's no significant improvement in likelihood by fitting separate models, implying that a single kinetic model is sufficient to explain the gene's dynamics across the tested groups.

# Why is this Powerful?

## Identifies "Kinetic Regimes"
It allows scVelo to pinpoint genes whose regulatory mechanisms change fundamentally between cell types or lineages, providing insight into distinct "kinetic regimes."

## Refined Trajectory Inference
By recognizing and accounting for these differential kinetics, scVelo can build more accurate and nuanced developmental trajectories, avoiding the distortion that would arise from forcing a single kinetic model where it doesn't apply.

## Discovery of Regulatory Differences
Genes flagged by the differential kinetic test are prime candidates for further investigation into underlying regulatory differences (e.g., cell-type specific transcription factor binding, alternative splicing/polyadenylation events, or differential degradation machinery).

## Higher Resolution Insights
It moves beyond simply asking "is this gene on or off?" to "how is this gene's dynamic behavior being regulated differently across cell states?"

# Mathematical Framework

The complete differential kinetics test can be formalized as:

{{< math >}} $$H_0: \theta_{\text{all clusters}} = \theta_0 = (\alpha_0, \beta_0, \gamma_0)$$ {{< /math >}}

{{< math >}} $$H_1: \theta_k = (\alpha_k, \beta_k, \gamma_k) \text{ for each cluster } k$$ {{< /math >}}

Under {{< math >}} $H_0$ {{< /math >}}, the test statistic follows:

{{< math >}} $$\text{LR} \sim \chi^2_{df}$$ {{< /math >}}

where {{< math >}} $df = 3 \times (K-1)$ {{< /math >}} and {{< math >}} $K$ {{< /math >}} is the number of clusters.

The p-value is then computed as:

{{< math >}} $$p\text{-value} = P(\chi^2_{df} > \text{LR})$$ {{< /math >}}

Genes with {{< math >}} $p\text{-value} < \alpha$ {{< /math >}} (typically {{< math >}} $\alpha = 0.05$ {{< /math >}}) are considered to have significantly different kinetics across clusters.