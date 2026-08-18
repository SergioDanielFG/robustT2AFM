# Design decisions

Decisions that cannot be inferred from reading the code, together with the
measurements that support them. Each section answers a question that has been
asked more than once, almost always from the wrong direction.

---

## 1. Masking makes the classical chart deaf, not jumpy

This is the easiest thing in the whole method to get backwards, and it is worth
writing down because the natural intuition points the wrong way.

The wrong intuition says: if Phase 1 is contaminated, the classical chart will
fire constantly. The actual chain runs the other way:

1. Contamination in Phase 1 inflates the pooled covariance `Sp`. On the
   Tennessee Eastman data its determinant is **16.14 times** that of the robust
   reference matrix `Sw`.
2. The classical control limit **does not move**. `hotelling_classical_ucl()`
   is parametric and depends only on K, I, J and alpha: it returns 19.4644 with
   a clean Phase 1 and with a contaminated one.
3. What changes is the statistic. An inflated `Sp` shrinks every batch T² at
   once, and shrunken statistics never reach a limit that has stayed put.

Contamination does not make the classical chart jumpy. It makes it deaf. What
any figure or passage should show is **lost detection**.

The acceptance ellipsoid does widen in data space, which is where the loose
intuition comes from. But the UCL is a specific number and that number does not
change. See also section 10.

### The measurements

On `afm_phase2`, the T² of `F2_B01` drops from 39.19 under clean calibration to
16.92 under contaminated calibration, and the classical alarm count falls from
9 of 20 to 2 of 20, against the same limit of 19.4644.

Across the Tennessee Eastman application and the package example scenario, the
classical method produced **zero false alarms** in every case measured:

- Table 9 of the paper.
- The five faults of Table 10.
- The 20 Phase 2 compositions of the stability analysis.
- The 12 seeds and 200 replicates of the package example scenario.

Not "few". Zero.

### Two necessary qualifications

**In simulation the classical method does raise false alarms.** Table 1 of the
paper reports 146 over 200,000 in-control batches, FAR = 0.000730,
ARL0 = 1370. The zero holds for the Tennessee Eastman application and for the
package scenario; it is not a general property of the classical chart.

**And the proposed method is not free of them either.** When Phase 1 is
calibrated on IDV 1 batches, a far more severe fault than IDV 7, the reference
center drifts by 0.82 standard deviations and six of the ten fault-free Phase 2
batches cross the limit. This is Section 4.6 of the paper. The weighting
protects `Sw`, not `mu_r`.

### What follows from this

No figure, help page or vignette should suggest that the classical method
raises more false alarms under contamination. And where the count is zero, that
zero is not a shortcoming of the chosen scenario to be fixed by finding a
different one. Going looking for such a scenario is the inverted intuition all
over again.

This is what `plot_method_comparison()` exists for, and why its caption opens
with the mechanism rather than with a count.


---

## 2. Colouring by truth and colouring by verdict

Figure 6 of the paper colours batches by **truth**: fault-free against faulty.
That choice is what makes the argument visible, and it is legitimate because the
Tennessee Eastman collection **comes labelled** — we know which batches carry a
fault because the dataset says so. The resulting image is the whole argument:
twenty red points, nearly all of them below the limit in the classical panel.
Batches known to be defective, and a chart that says nothing.

A control chart in production cannot colour this way, and not for a technical
reason: **which batches carry a fault is precisely what we are trying to find
out**. If that information were available, the chart would be redundant.

Hence the rule the package follows:

- `plot_control_chart()` always colours by **verdict** — above or below the
  limit — and never uses contamination labels. This is the production chart.
- `plot_method_comparison()` colours by **truth** only when the caller supplies
  that truth, through the `faulty` argument. Without it, it falls back to
  colouring by verdict. Asking for truth colouring without supplying the truth
  is an error, not a silent default.

Colouring by truth is an instrument of **validation**, not of operation. It
serves to show why one method detects and another does not, on data where the
answer is known in advance. It is not a way to monitor a process.

Historical note: for a while the package could not reproduce its own Figure 6,
because the only chart function coloured by verdict by design and the published
figure came out of an external script. `plot_method_comparison()` with `faulty`
closes that gap without relaxing the rule above.

---

## 3. The vocabulary: two pairs, and why there are two

The paper uses two pairs of terms, and it is right to do so, because they name
two different things: "out of control" is what **the chart says**, and "faulty
batch" is what **the batch is**. Verdict against truth, which is the same
distinction that underpins the design of the charts in section 2.

The problem was never having two pairs. It was having five ways of saying them.
The package fixes these:

| Idea | Pair | Legend title |
|---|---|---|
| Chart's verdict | **Out of control / In control** | `Chart verdict` |
| Batch's condition | **Faulty batch / Fault-free batch** | `Batch condition` |

**Why "in control" and not "under control".** *Under control* is a calque from
Spanish. The statistical process control literature says *in control* and *out
of control*, without exception in Montgomery. The package console already said
"out of control" and "in control" in `print` and `summary`, so this choice
required no rewriting.

**Why "fault-free" and not "normal".** *Batch under normal conditions*
translates naturally as *normal batch*, and in this particular paper that is an
expensive ambiguity: the text is full of normal distributions. *Fault-free* is
the exact complement of *faulty*, contains no form of the word *control*, and
describes exactly what happens in the Tennessee Eastman data, where a batch
either carries an injected fault or it does not.

**Why the two pairs can coexist.** Three conditions, and all three hold. They
never appear in the same legend, because the colour encodes one thing or the
other and never both. They share no discriminating word, only the noun *batch*.
And the legend title names which of the two is being shown. That third condition
is what really settles it: a reader who sees two figures in a row does not have
to wonder whether the same colour means the same thing, because the legend says
so.

### The third level: the `Status` column of the datasets

The same decision had to be applied one level down. `simulate_batch_process()`
was generating the `Status` column with the values "Under Control" and "Out of
Control", and that column is **ground truth**, not a verdict: it records which
batches were generated with a fault.

The symptom was the documentation example, which built the `faulty` argument
like this:

```r
faulty <- afm_phase2$Batch[afm_phase2$Status == "Out of Control"]
```

An argument named `faulty`, which is truth, filtered by a label that names a
verdict. The example was teaching exactly the confusion the figures had just
removed.

`Status` therefore became **`Faulty` / `Fault-free`**, the same pair as the truth
legend. This was done **before the data were published**, and that timing
matters: once the `.rda` files are on Zenodo or CRAN, changing the values of a
column breaks the code of anyone already using them.

Verified by regenerating from the same seed: **columns `Var1` to `Var4` are
byte-identical** — `identical()` on the serialised objects, maximum difference
zero — the batch identifier does not move, and the three anchors still come out
at 19.69285, 13 and 19.46440, with a maximum T² of 44.8193. This was a change of
labels, not of data.

The same pass settled `ContaminationType`, which carried the level `OOC` in
Phase 2 — the same verdict jargon, one level down. It became **`Shifted`**, the
name Phase 1 already used, leaving the levels as
`Clean | Outliers | Shifted`.

The reason is that **it is literally the same line of code**: Phase 1 does
`mu + shift_contam * sigma_vec` and Phase 2 does `mu + shift_ooc * sigma_vec`,
and both go on to the same `gen_batch(I, mu_k, Sigma)`. Same mean shift, same
covariance; only the magnitude differs, and magnitude is a parameter, not a
category.

There is a real difference between the two cases, but it is one of role rather
than of mechanism: in Phase 1 a shifted batch is contamination the method must
absorb, and in Phase 2 it is the signal the method must detect. That role is
determined entirely by the `Phase` column, so encoding it again in
`ContaminationType` was duplicating information — and duplicating it under two
different words invited the reader to believe they were two phenomena.

Each column now answers exactly one question: `ContaminationType` how the batch
was spoiled, `Phase` where it sits, `Status` whether it carries a fault.

### Carried over into the translation of the paper

Figure 6 is currently labelled "In-control batch" / "Faulty batch". In the
English manuscript **that legend must become "Fault-free batch" / "Faulty
batch"**, not the literal translation.

The reason: Figure 6 colours by **truth**, and "in-control batch" would use for
truth the same phrase the rest of the paper uses for the **verdict**. One phrase
would name two different things in adjacent figures, which is exactly the
confusion the split vocabulary avoids. It is a two-word change and an easy one
to lose in translation, which is why it is recorded here.

---

## 4. Why the AFM weight plot does not colour by the 1/K line

`plot_afm_weights()` draws the reference at the uniform weight `1/K` but does
not colour the bars according to whether they fall below or above it. The
temptation is obvious and the reason to resist it is arithmetic.

The AFM weights are inverse and normalised:


They sum to one by construction, so their mean is exactly `1/K`. And the
distribution of `1/lambda1` is right-skewed, so its median sits below its mean.
The consequence: **more than half the batches fall below `1/K` even with a
clean Phase 1**.

Measured on the package scenarios with K = 30: with 2 genuinely contaminated
batches, **17 of 30** fall below `1/K`; with 6 contaminated batches, also **17
of 30**. Colouring by that threshold would paint more than half of a sound
calibration as suspect, and would communicate something false with a great deal
of confidence.

The bars are therefore drawn in a neutral colour, and emphasis is left to the
caller through `highlight_lowest`, which is a presentation choice rather than a
judgement made by the package.

**And this is why its reference line is not the alarm colour.** In both control
charts the line is a **limit**: crossing it raises an alarm, and the dark red
`#A02D31` is the alarm colour throughout the package. Here `1/K` **is not a
limit**: crossing it means nothing, as the 17 of 30 just showed. It is drawn in
slate `#2C3E50` to say that it is a reference and not a threshold. The two
colours distinguish two kinds of line; it is not an oversight, and it is written
in the `@details` of all three functions.

It is worth restating what the weight measures: internal dispersion, through the
first eigenvalue. It does not measure position. A batch translated as a whole,
with its shape unchanged, keeps its weight intact. The weighting protects `Sw`;
it does not protect the reference center `mu_r`. In the base scenario, of the six
lowest weights four belong to batches that do carry outliers and two to clean
batches: the weight ranks dispersion, it does not classify batches.


---

## 5. Phase 1 batches of unequal size

The control limit assumes a common batch size I, through `m* = round(I·h)`. The
calibration itself does not need one: MCD, the weights, `Sw` and `mu_r` are all
computed batch by batch and handle different sizes without difficulty. The limit
does need one.

In plant data unequal batches are the norm — one run stops early, another loses
discarded measurements. When that happens there is no exact I, and one has to be
chosen.

**The size chosen is the one whose `m*` equals the mean of the `m*_k`.** With
unequal batches the actual degrees of freedom of the covariance estimator are
`sum_k (m*_k - 1)`; setting that equal to the `K(m* - 1)` of the published
formula gives `m* = mean(m*_k)`, and hence:


**This is not the same as rounding the mean of the sizes.** The two are easy to
confuse and they genuinely differ:

- `mean(round(I_k · h))` is the one that follows from the degrees of freedom.
- `round(mean(I_k) · h)` is the one that comes from averaging sizes first and
  converting afterwards.

Over 20,000 random vectors of 30 sizes between 10 and 25, the two disagree in
**18.6%** of cases. A worked example that can be checked by hand, with h = 0.67
and sizes 11, 11, 12, 20, 20, 20:


In the worst divergence found (30 batches between 10 and 25) `m*` shifts from 11
to 12 and the UCL from 19.81864 to 19.74987, a difference of 0.347%. With equal
batches the two forms always agree, so nothing already measured moves.

The code uses the first form, the derived one. The rounding is lossless in both
directions: `round(round(m/h) · h) = m` for every m between 5 and 30.

**The minimum is not used**, which is what instinct calls the conservative
choice. The UCL decreases monotonically with I: with K = 30, J = 4 and h = 0.67,
I = 14 gives 20.0097, I = 18 gives 19.7499, I = 20 gives 19.6929 and I = 25
gives 19.5374. The minimum gives the **widest** limit, and a wider limit signals
less often. That would be prudent against false alarms, which in the Tennessee
Eastman application are not the binding risk, at the cost of detection, which
is. See section 1.

The magnitude of the difference is small, around 1.6% between I = 14 and I = 20.
What mattered was not the size of the error but that it was arbitrary: the
earlier implementation took the size of the **first** batch, and because batches
are traversed in order of appearance, **reordering the rows of the data frame
changed the limit** without changing a single value. `unique(data$Batch)` takes
order of appearance, and there is now a test that pins this.

`calibrate_afm_mcd` returns `batch_sizes` for this reason; without it
`ucl_F_adjusted` cannot check the assumption it depends on.

**Verified on `afm_phase1`.** The robust center comes out at
0.053 / 0.014 / 0.036 / 0.011 with 6 of 30 batches contaminated. The classical
center from `hotelling_classical_calibrate` on the same data comes out at
0.167 / 0.142 / 0.157 / 0.140 — 4.6 times further from the true zero.

---

## 6. Why the example uses seed 20260425

All of the package documentation — examples, vignette and README — uses a single
scenario, the base configuration of the paper:

```r
simulate_batch_process(
  K1 = 30, K2 = 20, I = 20, J = 4, rho = 0.6,
  outlier_batches_F1 = 6, outlier_rate = 0.20, outlier_shift = 4,
  prop_ooc_F2 = 0.5, shift_ooc = 1.0,
  seed = 20260425
)
```

The previous scenario had a deeper flaw: both methods flagged the same 6 of 20
batches and the determinant ratio was 1.02. Anyone running the default example
could not see what the method is for.

**The seed is deliberate and it is not the most flattering one.** Over 200
replicates of this configuration the package averages 0.824 detection for the
proposed method and 0.208 for the classical one. Seed 20260425 gives 8 of 10 and
2 of 10, the single replicate closest to both means at once. Seeds 20260417 and
20260419 give the proposed method 10 of 10: they are the most favourable outcome
and are ruled out for that reason. The example should not be "improved" by
switching to them.

**The example counts are not rates, and they do not measure the same thing as
the paper.** The paper reports 0.863 and 0.270, averaged over 2000 repetitions,
and measures them after equalising both methods to a common ARL0 through the
empirical quantile of each method's own null distribution over 5000 in-control
batches. The example applies each method's operating limit instead. The two
quantities are close but not identical: the same ordering and the same
separation are to be expected, agreement to the decimal is not. The gap between
0.824 and 0.863, and between 0.208 and 0.270, is explained by that difference
but **has not been verified**, and no attempt should be made to close it by
adjusting the example.

With delta = 1.0 and alpha = 0.001, moreover, detection power is the only thing
separating the two methods in this scenario. Neither the 12 seeds nor the 200
replicates produced a single false alarm for either method, which is consistent
with section 1 and is not a limitation of the example.

### Further choices in `plot_afm_weights`

**Transposed axes.** Batch identifiers go on the y axis so that they stay
horizontal and legible with the K = 30 batches of the paper. Sorting by weight
already discards run order, so nothing is lost by transposing.

**Observed range of the weights.** On `afm_phase1`: from 0.0065 (F1_B17) to
0.0900 (F1_B04), a factor of nearly 14 between the extremes with only 6 lightly
contaminated batches.

### Design choices in `plot_method_comparison`

**Why the figure carries no caption.** The central message is that the classical
chart does not become jumpy under contamination, it becomes deaf. That sentence
is the entire point of the figure, and it is not printed underneath: a caption
nobody reads is worse than none, because the lines that do matter get skipped
along with it.

**Why `scale = "ratio"` is the default.** With raw statistics the two panels need
independent vertical axes, and then their heights are not comparable. That calls
for a warning, and this figure has no caption to put one in. Dividing by each
method's own limit removes the need for the warning rather than repeating it.
With `scale = "T2"` the warning line reappears automatically: the caption exists
only when the figure needs defending against its own natural reading.

**Verified on `afm_phase2`.** The same batch F2_B01 sits at 1.83 of the limit in
the robust panel and at 0.87 in the classical one. F2_B04 likewise: 1.79 against
0.87. The classical alarm set is contained in the robust one, without exception.
It is the demonstration of masking without a single formula.

### Names withdrawn from the documentation of `ucl_F_adjusted`

**"Way 1".** The help page used this name for the analytic limit, citing it as
though the paper used it. It does not: Section 2.5 says "analytic limit" and
"operating limit". The name came from an earlier stage in which the candidate
routes were numbered (analytic, K_eff, bootstrap). When the others were
discarded the name lost its meaning and survived as a fossil. Replaced by a
reference to Equation (8).

**Hardin-Rocke.** The help page stated that this correction is not applied.
Hardin and Rocke appear neither in the paper nor in its bibliography; the
mention came from a separate clarification document. A denial with no reference
and no explanation informs nobody: the reader learns neither what is not being
applied nor why it should matter. Withdrawn.

**The two alphas.** `mcd_alpha` (the MCD retention fraction, 0.67) and `alpha`
(the false alarm rate, 0.001) are unrelated parameters with similar names. The
signatures are published and do not change; the warning goes in the help text of
both arguments.

### The `sample()` call in `simulate_batch_process`

**The bug.** `sample(x, n)` draws from `1:x` when `x` has length one, instead of
from `x` itself. It is a historical convenience in R that turns into a defect as
soon as a vector is down to a single element. Here it happened with
`available_for_outliers`: with `prop_contam_F1 = 0.97` and `K1 = 30`, 29 batches
are shifted and one remains available, and the outliers were then injected into a
different batch from the one chosen, leaving the `Status` and
`ContaminationType` labels misassigned with no warning at all.

**The fix.** `x[sample.int(length(x), n)]`, the idiom recommended in `sample`'s
own help page.

**Verified not to move the datasets.** With `prop_contam_F1 = 0`, which is the
configuration of `afm_phase1`, no earlier branch consumes the generator and
`available_for_outliers` arrives with all 30 batches, so both forms consume
randomness identically. The 24 tests in `test-data.R`, which regenerate the data
from the seeded call and compare against what is shipped, all still pass.

**The seed in that example.** `seed = 20260417` happens to be one of the two
seeds ruled out for the paper scenario for giving 10 of 10. The scenario here is
a different one (K1 = 30 with 2 outlier batches and 7% shifted, K2 with 30% out
of control), so it is not the same case, but the number coincides and can confuse
anyone reading both.

---

## 7. What each invalid input actually does, measured

The validation guards were added against a list of assumed failures. Measuring
them before writing the code showed that three of the four behaved differently
from what had been assumed, and in two cases worse. They are recorded here
because intuition about them has been wrong, and because they explain why Phase
1 and Phase 2 need different guards rather than one shared guard.

### A text column does not make `covMcd` fail: the calibration lies

A non-numeric variable was assumed to make `covMcd` fail with an internal
message from `robustbase`. **It does not fail at all.** `data.matrix()` converts
the text column into factor codes, and the levels are ordered alphabetically,
not numerically. On `F1_B01`:


The calibration finishes without a warning and returns `mu_r["Var2"] = 10.32`
and `Sw["Var2","Var2"] = 44.09`, against roughly 0.01 and 1.4 in the honest
columns. All finite, all plausible, all wrong. It is the same failure mode that
motivates the variable auto-detection message: a credible and incorrect result,
which is worse than an error.

### The text is silent only in Phase 1, and the reason matters

In Phase 2 the same text column **does** abort, with `'x' must be numeric`. The
difference is which function touches the data: Phase 1 goes through `covMcd`,
which coerces via `data.matrix()`; Phase 2 goes through `colMeans()`, which
refuses.

The two phases therefore do not carry the same guard out of symmetry. They carry
different guards because they fail differently: in Phase 1 the type check
prevents a silent disaster, and in Phase 2 it only improves a message that was
already loud.

### An NA in Phase 2 does not hide one batch: it destroys the whole count

`T2` comes out NA, so `is_ooc` comes out NA, and because the count is
`sum(mon$is_ooc)`, **the entire sum comes out NA rather than a number**. It is
not one batch's alarm that is lost: it is the total. This happened identically in
both twins, robust and classical, and reached the user in print through the real
path, `run_afm_mcd(compare_classical = TRUE)` followed by `summary()`.

This is why the Phase 2 guard was applied to both twins, whereas the rejection
of J = 1 was applied only to the main one: J = 1 is an absurd input with a
visible result, while an NA appears only in production and with an invisible
effect.

### An NA in Phase 1 is absorbed by `covMcd`

This is the only one of the four that turned out **less** serious than assumed.
Injecting an NA into `afm_phase1`, `mu_r` comes out identical to four decimals
(0.0535 / 0.0137 / 0.0356 / 0.0113 either way): `covMcd` absorbs it without
propagating. The guard is kept as a precaution and for symmetry with the
classical twin, which already had one, not because it currently produces bad
numbers.


---

## 8. Known debts

Things that are wrong or incomplete on purpose, with the reason they were left
that way. These are not open questions awaiting investigation: they have been
measured and decided, and what remains is to do them.

None of them changes a number in the paper or affects the reproducibility of its
tables.

### `hotelling_classical_calibrate` accepts J = 1 and its twin does not

`calibrate_afm_mcd` rejects a single variable with a message pointing the user
towards a univariate Shewhart chart. The classical twin does not: it accepts
J = 1 and returns a degenerate result.

The asymmetry is deliberate for now. The classical function is the reference
implementation of Equation (1) of the paper, not a production tool, and its
degenerate result with J = 1 is visible rather than silent. It is still lopsided,
and the two should either be aligned or the reason documented.

### Two error messages say what failed but not what to do

They are `"The following 'variables' are not numeric: ..."` and `"Batch 'X'
contains non-finite values (NA/NaN/Inf)."`. They break the convention the rest of
the package follows, which is to say what to do and not only what failed.

They sit in all four calibration and monitoring functions because they were
copied verbatim from `hotelling_classical_calibrate` when the guards were added,
so that the same failure would produce the same message in all four places.

**Fix them in a single pass across all four, not one at a time.** The value of
having copied them verbatim is precisely that uniformity; a partial fix leaves
one twin better than the other and destroys the only thing that copying bought.

### Two classical functions have no dedicated tests

`test-hotelling_classical.R` covers `hotelling_classical_monitor` and the
delegation from `run_afm_mcd`. `hotelling_classical_calibrate` and
`hotelling_classical_ucl` have no tests of their own: they are exercised
indirectly from `test-batch_col.R` and `test-run_afm_mcd.R`, which is fine for
what those files are testing but is no substitute.

Of the two, the more exposed is `hotelling_classical_ucl`, because it produces a
published anchor (19.46440) and has the least coverage relative to what it
guarantees.

### The non-finite guard lives inside the loop on purpose

When the two `monitor_*` functions were vectorised — the T² values accumulate in
vectors and the data frame is assembled once, instead of an `rbind` per iteration
— an obvious temptation appeared: lift the finiteness check into a single pass
over the whole data frame, which is faster and looks equivalent.

**It is not equivalent, and the difference breaks nothing visible.** The error
names a batch, and with the guard inside the loop that batch is **the first
invalid one in order of appearance**. A prior pass over the whole data frame
would name the first by row, which is not the same batch when the rows do not
arrive grouped by batch. The result is still a correct error on invalid data;
only which of the bad batches gets named changes, and that is exactly what an
engineer uses to go and find the problem.

A test pins this by reordering the rows so that order of appearance changes, and
checking that the named batch changes with it. That is the test which catches
this optimisation if anyone makes it on instinct.

### The two twins treat a single-observation batch differently

`hotelling_classical_monitor` has a branch for `I == 1`
(`as.numeric(subset_batch[1, ])` instead of `colMeans()`), and `monitor_afm_mcd`
does not.

The asymmetry is real and pinned by a test, but it is neither explained in the
code nor justified in the paper. With a single observation the batch T² stops
being the statistic the method defines — the batch mean is the observation itself
and there is no averaging to reduce the variance — so the reasonable outcome is
that neither function should accept it silently, rather than one handling it
separately. Which of the two forms is correct remains to be decided; until then,
the test prevents the branch from disappearing through carelessness during
refactoring.

### `Rplots.pdf`
