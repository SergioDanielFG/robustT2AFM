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
