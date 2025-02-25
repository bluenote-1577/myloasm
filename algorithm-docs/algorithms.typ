#set heading(numbering: "1.")
#import "@preview/ctheorems:1.1.2": *
#show: thmrules.with(qed-symbol: $square$)

#set math.equation(numbering: "(1)")

#let proof = thmproof("proof", "Proof")
#let theorem = thmbox("theorem", "Theorem", fill: rgb("#eeffee"))
#let definition = thmbox("definition", "Definition", fill: rgb("#ffeeee"))
#let example = thmbox("example", "Example", fill: rgb("#fffeee"))
#align(center, text(17pt)[
  Algorithms for myloasm - living document
])

#align(center, text(17pt)[
  Jim Shaw - 2025
])

#set par(spacing: 1.0em)
#set par(leading: 0.8em, first-line-indent: 1.8em, justify: true)
#show heading: set block(above: 1.5em, below: 1.5em)

= Distances between unitig-read distributions

Let $X = (x_1,x_2,x_3,...)$ and $Y = (y_1,y_2,y_3,...)$ be two distributions of read depths for two unitigs. Notably, each $x_i in bb(R)^n$ are vectors given a range of SNPmer identity values ($S$).

Our goal is to calculate a distance between these two distributions: $D(X,Y)$. $D(X,Y)$. Given an edge, $e(X,Y)$, between a unitig path (represented by $X$) and an adjacent path or contig ($Y$), we use $D(X,Y)$ to tell how "favourable" this connection is. The "probability" of this connection is a function of the distributions $X,Y$ and also the overlap properties of the edge $e$. Therefore, 

$ Pr(e(X,Y)) = exp(-[D(X,Y) + f(e)] / T) $ given some temperature $0.5 <= T <= 2$. The probability of a path will be a product of edge probabilities.

= Defining $f(e)$. 

Given an overlap, we consider the overlap length $ell(e)$ and the overlap identity $sigma(e)$. For now, we let 

$ f(e) = [1 - ell(e) / (max_{e in X^+} {ell(e)})] * c(e) + [ max_{e in X^+} {sigma(e)} - sigma(e)] $

We will define $c(e)$ to be a coverage control factor for overlaps. Lower converages mean higher variance in overlap lengths, so the first term should be less confident at low coverage. For now, this is something like $1 - 3/(min(x[-1], y[-1]) + 3)$ where $x,y$ are the adjacent reads and $x[-1]$ is the last entry of the vector $x$, where more stringent overlap coverage thresholds have later entries.

The first term controls for the relative overlap length of the edge compared to other adjacent edges from $X$: $X^+$. We consider $e$ a directed edge. 

The second term controls for the relative difference in overlap identity compared to adjacent edges. The idea is that given 0.99 vs 1.0 $sigma$ values, we get a factor of $e^(-1/T)$. Similarly, a vastly shorter overlap gives a factor of $e^(-1/T)$. Given no adjacent edges, this score is of course 0. 

= Defining $D(X,Y)$ 

Given a vector $x$, define $ log_a (x[i]) = log(x[i] + a) $ as a log with pseudocount. We let $a = 3$ in general. Then define the coordinate-wise log distance as:

$ d(x^i,y^i) = abs(log_a (x[i]) - log_a (y[i])). $


Intuitively, if unitig $X$ has coverage 20 and unitig $Y$ has coverage 10, then given some scaling factor $C$ (that can depend on a range of factors), $D(X,Y)$ should be about $C log(20/10) = C log(2)$. This gives a probability of $e^(- C log(2) / T) = 1/2^(C/T)$. 


The above intuition should carry over to the _distribution_ over $X, Y$: a unitig will not just have "coverage 20". Our goal is to define a distance over the distribution $X,Y$ in a principled way so that 

1. unitigs from the same genome should have similar distributions, including variances and means, and

2. account for uncertainty in $X$ and $Y$, which are extremely noisy and have small sample sizes (e.g. a unitig with only one read) and high variance due to inter/intragenomic repeats

== Dealing with uncertainty and variance in $X$ and $Y$ 

Let $M(X)$ be the coordinate-wise median of the distribution $X$ and also similarly for $Y$. Given samples $X$ and $Y$ from their respective, unknown probability distributions, we try to reason with $d(M(X)^i, M(Y)^i)$. We also have the distribution of pairwise log differences, $d(X^i, Y^i)$ on hand. Let's assume $d(X^i, Y^i)$ has some average $mu$ and standard deviation $Sigma$. 

=== Variance matching

=== Variance

To estimate the standard deviation, we take the IQR of the pairwise log differences. Theoretically, many distributions have IQR $prop$ $Sigma$. 

=== Sample size

To deal with sample sizes $N_x$ for $X$ and  $N_y$ for $Y$, we proceed by normalizing by a scaled confidence interval length. Let CI be the confidence interval length for the median log ratio given some confidence %. Given a distribution $d(X^i, Y^i)$, we normalize by $1 + "CI"(d(X^i, Y^i))$ which should go to 1 as $N_x, N_y$ get large. 

We can calculate this confidence interval in many ways. We could do boostrapping, where we resample $d(X^i, Y^i)$ and take medians. The main issue is that when $Y$ consists of a single read. This is because the "sample size" is essentially 1 here, but the bootstrap "feels" like the sample size is $N_x$. 

Instead, we propose the following. Under the assumption of normality, the confidence interval length is $prop Sigma / sqrt(N)$. We take $N$ as $2/(1/N_x + 1/N_y) = H(N_x,N_y)$, the harmonic mean, to bias for lower sample sizes. We already have a robust estimator proportional to $Sigma$: the IQR. 

=== Final formulae

==== Version 1 (Didn't work out well)

$ D(X,Y) = sum_(i=1)^n d(M(X)^i, M(Y)^i) / (([0.5 + 1/(1+max(N_x,N_y))]+ hat(Sigma)) * (1 + hat(Sigma)/sqrt(H(N_x,N_y)))) $

Update: Didn't work well. Turns out we _want_ to penalize high variance: we want our final paths to have low variance. 

==== Version 2 (Worked okay; 2-20 tests "POST-IQR" code in notebook)

$ V_i = abs(d(U(X)^i, U(Y)^i) - d(L(X)^i, L(Y)^i)) $
$ D_i = d(M(X)^i, M(Y)^i) $

$ D(X,Y) = sum_(i=1)^n (D_i + V_i) / (1 + hat(Sigma)/(4 sqrt(H(N_x,N_y)))) $

Update: Worked okay. However, the sum formulation has issues when there is an interspecies repeat. First coordinate will have large variance but the more strict coverages will have low (i.e. good) variance. 

==== Version 3 (TODO)

$ D_i = d(M(X)^i, M(Y)^i) $
$ V_i = abs(d(U(X)^i, U(Y)^i) - d(L(X)^i, L(Y)^i)) $

$ D(X,Y) = min_(i) (D_i + V_i) /  (1 + hat(Sigma)/(4sqrt(H(N_x,N_y)))) $





