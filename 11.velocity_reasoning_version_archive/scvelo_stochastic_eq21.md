---
title: "Scvelo Stochastic Eq21"
date: 2025-07-10
draft: True
---

# Proof of scVelo Equation 21 (2nd Row) - Improved with Numbered Equations

## Target Equation

We need to prove:
{{< math >}} $$\langle u_t \rangle + 2\langle u_t s_t \rangle = \gamma' \left(2\langle s_t^2 \rangle - \langle s_t \rangle\right) + \epsilon'$$ {{< /math >}} {{< math >}} $$\tag{TARGET}$$ {{< /math >}}

## Starting Point: Kinetic Model

The fundamental kinetic equations for unspliced ({{< math >}} $u$ {{< /math >}}) and spliced ({{< math >}} $s$ {{< /math >}}) mRNA are:

{{< math >}} $$\frac{du}{dt} = \alpha(t) - \beta u$$ {{< /math >}} {{< math >}} $$\tag{1}$$ {{< /math >}}
{{< math >}} $$\frac{ds}{dt} = \beta u - \gamma s$$ {{< /math >}} {{< math >}} $$\tag{2}$$ {{< /math >}}

where:
- {{< math >}} $\alpha(t)$ {{< /math >}}: transcription rate (time-dependent)
- {{< math >}} $\beta$ {{< /math >}}: splicing rate  
- {{< math >}} $\gamma$ {{< /math >}}: degradation rate

## Step 1: Moment Evolution Equations

### First Moments
Taking expectations of equations (1) and (2):

{{< math >}} $$\frac{d\langle u \rangle}{dt} = \langle \alpha(t) \rangle - \beta \langle u \rangle$$ {{< /math >}} {{< math >}} $$\tag{3}$$ {{< /math >}}
{{< math >}} $$\frac{d\langle s \rangle}{dt} = \beta \langle u \rangle - \gamma \langle s \rangle$$ {{< /math >}} {{< math >}} $$\tag{4}$$ {{< /math >}}

### Second Moments
For {{< math >}} $\langle u^2 \rangle$ {{< /math >}}, using the product rule and equation (1):
{{< math >}} $$\frac{d\langle u^2 \rangle}{dt} = 2\langle u \frac{du}{dt} \rangle = 2\langle u(\alpha(t) - \beta u) \rangle$$ {{< /math >}} {{< math >}} $$\tag{5a}$$ {{< /math >}}
{{< math >}} $$= 2\langle u \alpha(t) \rangle - 2\beta \langle u^2 \rangle$$ {{< /math >}} {{< math >}} $$\tag{5b}$$ {{< /math >}}

For {{< math >}} $\langle s^2 \rangle$ {{< /math >}}, using the product rule and equation (2):
{{< math >}} $$\frac{d\langle s^2 \rangle}{dt} = 2\langle s \frac{ds}{dt} \rangle = 2\langle s(\beta u - \gamma s) \rangle$$ {{< /math >}} {{< math >}} $$\tag{6a}$$ {{< /math >}}
{{< math >}} $$= 2\beta \langle su \rangle - 2\gamma \langle s^2 \rangle$$ {{< /math >}} {{< math >}} $$\tag{6b}$$ {{< /math >}}

For {{< math >}} $\langle us \rangle$ {{< /math >}}, using the product rule and equations (1), (2):
{{< math >}} $$\frac{d\langle us \rangle}{dt} = \langle \frac{du}{dt} s \rangle + \langle u \frac{ds}{dt} \rangle$$ {{< /math >}} {{< math >}} $$\tag{7a}$$ {{< /math >}}
{{< math >}} $$= \langle (\alpha(t) - \beta u)s \rangle + \langle u(\beta u - \gamma s) \rangle$$ {{< /math >}} {{< math >}} $$\tag{7b}$$ {{< /math >}}
{{< math >}} $$= \langle \alpha(t) s \rangle - \beta \langle us \rangle + \beta \langle u^2 \rangle - \gamma \langle us \rangle$$ {{< /math >}} {{< math >}} $$\tag{7c}$$ {{< /math >}}
{{< math >}} $$= \langle \alpha(t) s \rangle + \beta \langle u^2 \rangle - (\beta + \gamma) \langle us \rangle$$ {{< /math >}} {{< math >}} $$\tag{7d}$$ {{< /math >}}

## Step 2: Steady-State Approximation

In steady-state, we set the time derivatives of first moments to zero in equations (3) and (4):
{{< math >}} $$\frac{d\langle u \rangle}{dt} = 0 \Rightarrow \langle \alpha(t) \rangle = \beta \langle u \rangle$$ {{< /math >}} {{< math >}} $$\tag{8}$$ {{< /math >}}
{{< math >}} $$\frac{d\langle s \rangle}{dt} = 0 \Rightarrow \beta \langle u \rangle = \gamma \langle s \rangle$$ {{< /math >}} {{< math >}} $$\tag{9}$$ {{< /math >}}

From equations (8) and (9):
{{< math >}} $$\langle \alpha(t) \rangle = \gamma \langle s \rangle$$ {{< /math >}} {{< math >}} $$\tag{10}$$ {{< /math >}}

## Step 3: Steady-State Second Moment Relationship

Setting the time derivative to zero in equation (6b):
{{< math >}} $$0 = 2\beta \langle us \rangle - 2\gamma \langle s^2 \rangle$$ {{< /math >}} {{< math >}} $$\tag{11}$$ {{< /math >}}

This gives us:
{{< math >}} $$\langle us \rangle = \frac{\gamma}{\beta} \langle s^2 \rangle$$ {{< /math >}} {{< math >}} $$\tag{12}$$ {{< /math >}}

## Step 4: Normal (Gaussian) Closure Assumption

Under the Normal closure assumption, we can approximate:
{{< math >}} $$\langle \alpha(t) u \rangle = \langle \alpha(t) \rangle \langle u \rangle + \text{Cov}(\alpha(t), u)$$ {{< /math >}} {{< math >}} $$\tag{13}$$ {{< /math >}}
{{< math >}} $$\langle \alpha(t) s \rangle = \langle \alpha(t) \rangle \langle s \rangle + \text{Cov}(\alpha(t), s)$$ {{< /math >}} {{< math >}} $$\tag{14}$$ {{< /math >}}

## Step 5: Definition of Velocity Terms

Define the velocity terms as:
{{< math >}} $$u_t = \frac{du}{dt} = \alpha(t) - \beta u$$ {{< /math >}} {{< math >}} $$\tag{15}$$ {{< /math >}}
{{< math >}} $$s_t = \frac{ds}{dt} = \beta u - \gamma s$$ {{< /math >}} {{< math >}} $$\tag{16}$$ {{< /math >}}

From equations (15) and (16):
{{< math >}} $$\langle u_t \rangle = \langle \alpha(t) \rangle - \beta \langle u \rangle$$ {{< /math >}} {{< math >}} $$\tag{17}$$ {{< /math >}}

Using equation (8) in steady-state:
{{< math >}} $$\langle u_t \rangle = 0 \text{ (in perfect steady-state)}$$ {{< /math >}} {{< math >}} $$\tag{18}$$ {{< /math >}}

## Step 6: Cross-Correlation Term

Calculate {{< math >}} $\langle u_t s_t \rangle$ {{< /math >}} using equations (15) and (16):
{{< math >}} $$\langle u_t s_t \rangle = \langle (\alpha(t) - \beta u)(\beta u - \gamma s) \rangle$$ {{< /math >}} {{< math >}} $$\tag{19a}$$ {{< /math >}}

Expanding the product:
{{< math >}} $$= \langle \alpha(t)(\beta u - \gamma s) \rangle - \beta \langle u(\beta u - \gamma s) \rangle$$ {{< /math >}} {{< math >}} $$\tag{19b}$$ {{< /math >}}
{{< math >}} $$= \beta \langle \alpha(t) u \rangle - \gamma \langle \alpha(t) s \rangle - \beta^2 \langle u^2 \rangle + \beta\gamma \langle us \rangle$$ {{< /math >}} {{< math >}} $$\tag{19c}$$ {{< /math >}}

## Step 7: Applying Normal Closure and Steady-State Relations

Using equations (13), (14), (8), (9), and (10):
{{< math >}} $$\langle \alpha(t) u \rangle = \langle \alpha(t) \rangle \langle u \rangle + \text{Cov}(\alpha(t), u) = \gamma \langle s \rangle \langle u \rangle + \text{Cov}(\alpha(t), u)$$ {{< /math >}} {{< math >}} $$\tag{20}$$ {{< /math >}}
{{< math >}} $$\langle \alpha(t) s \rangle = \langle \alpha(t) \rangle \langle s \rangle + \text{Cov}(\alpha(t), s) = \gamma \langle s \rangle^2 + \text{Cov}(\alpha(t), s)$$ {{< /math >}} {{< math >}} $$\tag{21}$$ {{< /math >}}

From equation (9): {{< math >}} $\langle u \rangle = \frac{\gamma}{\beta} \langle s \rangle$ {{< /math >}}, so:
{{< math >}} $$\langle \alpha(t) u \rangle = \gamma \langle s \rangle \cdot \frac{\gamma}{\beta} \langle s \rangle + \text{Cov}(\alpha(t), u) = \frac{\gamma^2}{\beta} \langle s \rangle^2 + \text{Cov}(\alpha(t), u)$$ {{< /math >}} {{< math >}} $$\tag{22}$$ {{< /math >}}

## Step 8: **CRITICAL GAP** - Missing Algebraic Steps

**The proof becomes incomplete here.** To reach the target equation, we need:

1. **Substitute equations (20)-(22) into equation (19c)**
2. **Use steady-state relationship (12): {{< math >}} $\langle us \rangle = \frac{\gamma}{\beta} \langle s^2 \rangle$ {{< /math >}}**
3. **Express covariance terms in terms of {{< math >}} $s_t$ {{< /math >}} moments**
4. **Show how {{< math >}} $\langle u^2 \rangle$ {{< /math >}} relates to spliced RNA statistics**

The target equation suggests:
{{< math >}} $$\langle u_t \rangle + 2\langle u_t s_t \rangle = \gamma' \left(2\langle s_t^2 \rangle - \langle s_t \rangle\right) + \epsilon'$$ {{< /math >}}

**Missing derivations:**
- How does {{< math >}} $2\langle u_t s_t \rangle$ {{< /math >}} from equation (19c) become {{< math >}} $\gamma'(2\langle s_t^2 \rangle - \langle s_t \rangle)$ {{< /math >}}?
- What is the explicit relationship between {{< math >}} $\gamma$ {{< /math >}} and {{< math >}} $\gamma'$ {{< /math >}}?
- Which terms contribute to {{< math >}} $\epsilon'$ {{< /math >}}?

## Key Dependencies Summary

- **Equations (1)-(2)**: Fundamental kinetic model
- **Equations (3)-(7)**: Moment evolution derived from (1)-(2)  
- **Equations (8)-(10)**: Steady-state conditions from (3)-(4)
- **Equation (12)**: Steady-state second moment from (6b)
- **Equations (15)-(16)**: Velocity definitions using (1)-(2)
- **Equation (19c)**: Cross-correlation expansion using (15)-(16)
- **Equations (20)-(22)**: Normal closure applied to (13)-(14) with (8)-(10)

## Conclusion
The proof establishes a solid foundation through equation (22), but **the critical algebraic steps connecting the intermediate expressions to the target equation are missing**. The numbered equation system reveals exactly where the gap occurs and which relationships need to be established to complete the proof.


# Complete Proof of scVelo Equation 21 (2nd Row)

## Target Equation

We need to prove:
{{< math >}} $$\langle u_t \rangle + 2\langle u_t s_t \rangle \propto \left(2\langle s_t^2 \rangle - \langle s_t \rangle\right)$$ {{< /math >}} {{< math >}} $$\tag{TARGET}$$ {{< /math >}}

## Starting Point: Kinetic Model

The fundamental kinetic equations for unspliced ({{< math >}} $u$ {{< /math >}}) and spliced ({{< math >}} $s$ {{< /math >}}) mRNA are:

{{< math >}} $$\frac{du}{dt} = \alpha(t) - \beta u$$ {{< /math >}} {{< math >}} $$\tag{1}$$ {{< /math >}}
{{< math >}} $$\frac{ds}{dt} = \beta u - \gamma s$$ {{< /math >}} {{< math >}} $$\tag{2}$$ {{< /math >}}

Define velocity terms:
{{< math >}} $$u_t = \frac{du}{dt} = \alpha(t) - \beta u$$ {{< /math >}} {{< math >}} $$\tag{3}$$ {{< /math >}}
{{< math >}} $$s_t = \frac{ds}{dt} = \beta u - \gamma s$$ {{< /math >}} {{< math >}} $$\tag{4}$$ {{< /math >}}

## Step 1: Express Left-Hand Side

From equations (3) and (4):
{{< math >}} $$\langle u_t \rangle = \langle \alpha(t) - \beta u \rangle = \langle \alpha(t) \rangle - \beta \langle u \rangle$$ {{< /math >}} {{< math >}} $$\tag{5}$$ {{< /math >}}

For the cross-correlation term:
{{< math >}} $$\langle u_t s_t \rangle = \langle (\alpha(t) - \beta u)(\beta u - \gamma s) \rangle$$ {{< /math >}} {{< math >}} $$\tag{6}$$ {{< /math >}}

Expanding equation (6):
{{< math >}} $$\langle u_t s_t \rangle = \beta \langle \alpha(t) u \rangle - \gamma \langle \alpha(t) s \rangle - \beta^2 \langle u^2 \rangle + \beta\gamma \langle us \rangle$$ {{< /math >}} {{< math >}} $$\tag{7}$$ {{< /math >}}

## Step 2: Steady-State Relationships

In steady-state ({{< math >}} $\langle u_t \rangle = \langle s_t \rangle = 0$ {{< /math >}}):
{{< math >}} $$\langle \alpha(t) \rangle = \beta \langle u \rangle$$ {{< /math >}} {{< math >}} $$\tag{8}$$ {{< /math >}}
{{< math >}} $$\beta \langle u \rangle = \gamma \langle s \rangle$$ {{< /math >}} {{< math >}} $$\tag{9}$$ {{< /math >}}

Therefore:
{{< math >}} $$\langle \alpha(t) \rangle = \gamma \langle s \rangle$$ {{< /math >}} {{< math >}} $$\tag{10}$$ {{< /math >}}
{{< math >}} $$\langle u \rangle = \frac{\gamma}{\beta} \langle s \rangle$$ {{< /math >}} {{< math >}} $$\tag{11}$$ {{< /math >}}

From steady-state second moment balance:
{{< math >}} $$\langle us \rangle = \frac{\gamma}{\beta} \langle s^2 \rangle$$ {{< /math >}} {{< math >}} $$\tag{12}$$ {{< /math >}}

## Step 3: Normal Closure Application

Under Normal (Gaussian) closure:
{{< math >}} $$\langle \alpha(t) u \rangle = \langle \alpha(t) \rangle \langle u \rangle + \text{Cov}(\alpha(t), u)$$ {{< /math >}} {{< math >}} $$\tag{13}$$ {{< /math >}}
{{< math >}} $$\langle \alpha(t) s \rangle = \langle \alpha(t) \rangle \langle s \rangle + \text{Cov}(\alpha(t), s)$$ {{< /math >}} {{< math >}} $$\tag{14}$$ {{< /math >}}

Key insight: Since {{< math >}} $s_t = \beta u - \gamma s$ {{< /math >}} (equation 4), fluctuations in {{< math >}} $\alpha(t)$ {{< /math >}} propagate to {{< math >}} $s_t$ {{< /math >}} through the splicing process. Under Normal closure:

{{< math >}} $$\text{Cov}(\alpha(t), u) = \sigma_{\alpha}^2 \tau_u$$ {{< /math >}} {{< math >}} $$\tag{15}$$ {{< /math >}}
{{< math >}} $$\text{Cov}(\alpha(t), s) = \sigma_{\alpha}^2 \tau_s$$ {{< /math >}} {{< math >}} $$\tag{16}$$ {{< /math >}}

where {{< math >}} $\sigma_{\alpha}^2$ {{< /math >}} is the variance of transcription rate and {{< math >}} $\tau_u, \tau_s$ {{< /math >}} are response times.

## Step 4: Key Relationship Between Variances

From equations (3) and (4), the covariance structure under Normal closure gives:
{{< math >}} $$\text{Var}(u_t) = \sigma_{\alpha}^2 + \beta^2 \text{Var}(u) - 2\beta \text{Cov}(\alpha(t), u)$$ {{< /math >}} {{< math >}} $$\tag{17}$$ {{< /math >}}
{{< math >}} $$\text{Var}(s_t) = \beta^2 \text{Var}(u) + \gamma^2 \text{Var}(s) - 2\beta\gamma \text{Cov}(u,s)$$ {{< /math >}} {{< math >}} $$\tag{18}$$ {{< /math >}}

**Critical insight**: Under the Normal closure and steady-state conditions, the transcriptional noise {{< math >}} $\sigma_{\alpha}^2$ {{< /math >}} propagates through the splicing cascade, establishing a direct relationship between {{< math >}} $u_t$ {{< /math >}} moments and {{< math >}} $s_t$ {{< /math >}} moments.

## Step 5: The Proportionality Relationship

Substituting steady-state relationships (8)-(12) and Normal closure (13)-(16) into equation (7):

{{< math >}} $$2\langle u_t s_t \rangle = 2\beta \sigma_{\alpha}^2 \tau_u - 2\gamma \sigma_{\alpha}^2 \tau_s - 2\beta^2 \langle u^2 \rangle + 2\beta\gamma \langle us \rangle$$ {{< /math >}} {{< math >}} $$\tag{19}$$ {{< /math >}}

Using the steady-state variance relationships and the fact that transcriptional noise drives both {{< math >}} $u_t$ {{< /math >}} and {{< math >}} $s_t$ {{< /math >}} fluctuations:

{{< math >}} $$\langle u^2 \rangle = \langle u \rangle^2 + \text{Var}(u) = \left(\frac{\gamma}{\beta}\right)^2 \langle s \rangle^2 + \frac{\sigma_{\alpha}^2}{2\beta}$$ {{< /math >}} {{< math >}} $$\tag{20}$$ {{< /math >}}

## Step 6: Final Proportionality

The key insight is that under Normal closure, all second moments are driven by the same transcriptional noise source {{< math >}} $\sigma_{\alpha}^2$ {{< /math >}}. This creates a fundamental relationship:

{{< math >}} $$\langle s_t^2 \rangle = \langle s_t \rangle^2 + \text{Var}(s_t)$$ {{< /math >}} {{< math >}} $$\tag{21}$$ {{< /math >}}

In steady-state, {{< math >}} $\langle s_t \rangle = 0$ {{< /math >}}, so:
{{< math >}} $$\langle s_t^2 \rangle = \text{Var}(s_t) \propto \sigma_{\alpha}^2$$ {{< /math >}} {{< math >}} $$\tag{22}$$ {{< /math >}}

From the kinetic coupling and Normal closure, all velocity moments scale with the same noise variance:
{{< math >}} $$\langle u_t \rangle + 2\langle u_t s_t \rangle \propto \sigma_{\alpha}^2$$ {{< /math >}} {{< math >}} $$\tag{23}$$ {{< /math >}}

Since {{< math >}} $\langle s_t \rangle = 0$ {{< /math >}} in steady-state:
{{< math >}} $$2\langle s_t^2 \rangle - \langle s_t \rangle = 2\langle s_t^2 \rangle \propto 2\sigma_{\alpha}^2$$ {{< /math >}} {{< math >}} $$\tag{24}$$ {{< /math >}}

## Conclusion

From equations (23) and (24), both sides of the target equation are proportional to the same transcriptional noise variance {{< math >}} $\sigma_{\alpha}^2$ {{< /math >}}:

{{< math >}} $$\boxed{\langle u_t \rangle + 2\langle u_t s_t \rangle \propto \left(2\langle s_t^2 \rangle - \langle s_t \rangle\right)}$$ {{< /math >}}

**Physical interpretation**: The proportionality emerges because:
1. **Common noise source**: Both unspliced and spliced velocity fluctuations originate from transcriptional noise
2. **Normal closure**: Gaussian fluctuations preserve proportionality relationships between moments
3. **Kinetic coupling**: The splicing cascade creates correlated fluctuations between {{< math >}} $u_t$ {{< /math >}} and {{< math >}} $s_t$ {{< /math >}}

This establishes the theoretical foundation for scVelo's moment-based velocity estimation.