# Node Census

## 1 What a Node Is

A **node** is a point where a quantity changes meaning, not where a variable changes value. A rename is not a node. Neither is a reshape.

A node is one of the following:

* An aggregation, join, or filter that changes what one row represents
* A call into a statistical package
* Arithmetic combining two derived quantities
* A filter that changes which population is described

Three rules settle the cases where the same code appears more than once.

1. **Same operation, same inputs, same settings:**
   Two steps performing the same operation, on the same inputs, with the same settings, are **one node**. The duplication is itself a finding, because nothing keeps two copies in agreement and a correction made at the node reaches only one of them.

2. **Same code, different inputs:**
   The same code applied to different inputs is **two nodes**, because the inputs are what is being checked.

3. **Function definition and call sites:**
   A function definition called from several places is **one node**, and each call site is an additional node.

   In our pipeline, one coverage-standardization function is called on three different assemblages. The function is one node, its three call sites are three more, and the difference between them is the substance of the analysis.

---

## 2 The Region a Node Occupies

Every node occupies a **contiguous region of the file**, delimited by a first and a last line copied word for word from the source.

The region, not the single line, is what a reviewer reads when the node comes up for checking.

A region will usually contain statements that are not themselves nodes, such as:

* A reshape
* A rename
* An intermediate object that carries a result to the next line

These are absorbed into the node whose decision they serve. The region is what makes the absorption visible.

Recording in one clause why a statement was absorbed keeps absorption from becoming a way of leaving code unread.

### Region Rules

* Regions must not overlap.
* Rows are listed in source order.
* A statement that falls outside every region is a claim that it belongs to no node, and that claim should be one you intend to make.

### Anchors, Not Line Numbers

A region is located by quoting its **first and last lines**, because line numbers drift with every edit and a line of code does not.

The rules that make an anchor survive editing and be findable by a script are:

* **Copy the line exactly from the source.**
  Include its indentation and any trailing comma or pipe. Do not retype it, tidy the spacing, rewrite `%>%` as `|>`, or add/remove a `dplyr::` prefix.

* **Do not include line numbers, file paths, or comment text.**
  A line whose only distinguishing feature is its comment cannot serve as an anchor.

* **Make each anchor unique within its file.**
  If it is not unique, pair it with a context line above it that makes the pair unique.

* **Prefer stable lines.**
  An assignment with a named function call is more stable than a long line that a formatter might wrap.

* **Avoid double quotes where possible.**
  A CSV has to escape them, and a hand-edited table is easy to break.

The last line of a region is often a closing line such as `)`, which is almost never unique. That is expected. Keep the true region and supply a context line rather than shortening the region to obtain a tidier ending.

---

## 3 The Node Table

Use **one row per node**, in a plain-text table kept under version control alongside the code.

### Columns

| Column          | Description                                                                                                                                     |
| --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`            | A stable name, `script.quantity`. Other rows refer to it, so it must not change once cited.                                                     |
| `script`        | The filename.                                                                                                                                   |
| `section`       | The existing section header in the script that the node falls under.                                                                            |
| `start_anchor`  | The first line of the node region, copied exactly from the source.                                                                              |
| `start_context` | Context line used when the start anchor is not unique.                                                                                          |
| `end_anchor`    | The last line of the node region, copied exactly from the source.                                                                               |
| `end_context`   | Context line used when the end anchor is not unique.                                                                                            |
| `inputs`        | The node IDs this node consumes, or the source file for a node reading a fixed input. This is what makes the table a graph.                     |
| `grain_in`      | What one input row means; e.g., gear × bay × season × year × species.                                                                           |
| `grain_out`     | What one output row means. Where `grain_in` and `grain_out` differ, the difference is the thing to check.                                       |
| `units_in`      | Input units, such as counts, individuals, Hill numbers, nats, or proportions.                                                                   |
| `units_out`     | Output units.                                                                                                                                   |
| `operation`     | One sentence describing what the node does and why.                                                                                             |
| `settings`      | Decisions embedded in the node, including thresholds, target levels, bounds, tie-breaks, seeds, and any argument that changes what is computed. |
| `na_behavior`   | What happens when the operation cannot produce a value, and whether the causes are distinguishable downstream.                                  |
| `checked`       | The date and method of checking.                                                                                                                |
| `verdict`       | One of `ok`, `finding`, or `question`.                                                                                                          |
| `notes`         | The text of any finding or question.                                                                                                            |

The `settings` column is particularly important because it can become **Methods text**.

The `na_behavior` column is required for every node. Where a function returns a missing value for several different reasons that arrive downstream as the same `NA`, say so explicitly because undistinguished absence is otherwise invisible.

The final three columns—`checked`, `verdict`, and `notes`—are left empty during enumeration.

---

## 4 Building the Table

### Pass One: Enumerate

Enumerate **one script at a time, in source order**.

Read the script, mark the regions, and fill every column except the last three.

Three rules govern what goes into the table:

1. **Describe what the code does, not what it should do.**

   Where the code and its comments disagree, record the disagreement in `operation` rather than choosing between them.

2. **Document absorbed statements.**

   Where a region absorbs a statement that is not itself a node, say in one clause why it was absorbed.

3. **Do not guess.**

   Where the meaning of a quantity cannot be determined from the code alone, write the question in `operation` rather than guessing.

   An unanswered question in the table is worth more than a plausible answer.

### Pass Two: Link, Then Prune Backwards

The `inputs` column turns the rows into a **graph**.

Start from each reported number and walk backward through its inputs until you reach a fixed input file.

This bounds the census:

> Everything a reported number depends on is enumerated, and anything the walk does not reach is a candidate for deletion rather than for checking.

Enumerating forwards through a file and pruning backwards from the results is deliberate.

Reading a script in its own order is what a person can actually do without holding the whole pipeline in mind. The backward walk keeps the census bounded by the manuscript rather than by the entire codebase.

### Check Anchors Mechanically

Check the anchors mechanically, not by eye.

A short script should confirm that:

* Every start and end anchor occurs in its file.
* Each anchor is unique or has a context line.
* The end anchor occurs at or after the start anchor.
* No two regions overlap.

Report any anchor that could not be made unique rather than leaving it to be discovered later.

### Delegation and Its Limit

Pass one is largely mechanical and can be delegated, including to a language model given the column definitions and rules above.

The prompt used for this is provided as **Supplementary Material [Y] (`node_table_prompt.md`)**. It states the rules of Sections 4.1–4.3 in the form a model needs, fixes the output format, and requires the model to report any anchor it could not make unique rather than adjusting a region to make the check pass.

Pass two and the checking in Section 4.5 are **judgment-based and cannot be delegated**.

The reason is the same one that governs hand computation below:

> If a single source produces both the code and the account it is checked against, any misunderstanding enters both sides and the check passes.

Therefore, a table produced this way is read against the code row by row before any node is checked.

A description that is wrong in the same way the code is wrong is worse than no description because it makes the error look confirmed.

---

## 5 Checking the Table: Forwards

Order the completed table so that **every node appears after its inputs**, and check in that order.

This matters more than it sounds.

Checked backwards, every node raises questions that are answered upstream, and the reviewer must hold those questions open while descending.

Checked forwards, each node's inputs are already understood when it is reached, and the only new thing is the transformation itself.

In our experience, the same nodes take substantially longer backwards than forwards. The forward pass also matches how errors propagate:

> **Verify the source → verify the transformation → verify the consumer.**

### Four Questions at Each Node

At each node, ask:

1. Does the code in the region do what `operation` says?
2. Is `grain_out` correct?

   * Read what the operation groups or joins by.
   * Confirm that one output row means what the column claims.
3. Are the `settings` the intended ones?

   * A setting that cannot be justified once written down is a finding.
4. Does `na_behavior` describe what actually happens, and is it acceptable?

### Independent Arithmetic Checks

Where the operation is arithmetic that can be evaluated on paper, do so.

Choose cases where the computation reduces to something independently checkable. For example:

* Two assemblages sharing no species drive a multiplicative partition to its ceiling.
* A rarefaction target at or above the observed total makes the draw a no-op and reduces the computation to the underlying index on a known incidence matrix.

Write the expected value down **before running anything**.

The expected value must come from outside the code being tested.

If the same source produces both the code and the expected answer, any misunderstanding enters both sides and the test passes.

Where AI assistance was used to write the analysis, it should not also be used to generate the values the analysis is checked against.

---

## 6 What the Census Makes Visible

A node table is a description of the analysis in one place, at a uniform level of detail, with every quantity's:

* Grain
* Units
* Settings

recorded alongside the others.

Assembling it has an effect that checking the code does not:

> Choices that were made separately, months apart, in different scripts, end up on adjacent rows where they can be compared.

### Example 1: A Divisor That Differs Between Two Levels

The multiplicative partition divides an assemblage diversity by the number of assemblages.

At community level, the divisor is the number of bays sampled.

At guild level, it is the number of bays where the guild was observed.

Each is defensible in isolation.

The census put them on two rows with the same operation text and different settings. This made the difference a **question rather than an accident of two authors on two days**.

### Example 2: A Rescaling Whose Denominator Is Itself a Result

Guild beta is rescaled by its ceiling so that guilds occupying different numbers of bays are comparable.

The ceiling is observed occupancy, which is itself changing over the study period.

Therefore, a guild that contracts its range alters its own denominator.

Recording that in one field, next to the note that the rescaled value is undefined for a guild in a single bay, turned a line of arithmetic into a methodological question about what the trend measures.

Neither question is answered by a verdict in the table.

What the census does is produce these questions in a form specific enough to argue about and record the choice actually made so that a reader can disagree with it knowingly.

---

## 7 What the Census Does Not Cover

A node census subsumes most of what would otherwise be separate activities.

Reading the computational core is what checking a node is, organized by **quantity rather than by file**.

Shared conventions become node attributes, so a category level used inconsistently appears as two nodes disagreeing about the same field. This is stronger than a text search.

Grain mismatch becomes structural rather than something a reviewer must remember to look for.

Two things remain outside the census. Both are inexpensive and should be kept alongside it.

### 1. Text Search of the Whole Codebase

The census is organized by quantity, and each quantity appears once.

If the same wrong value was typed into four scripts, one node describes it and the table shows one instance. Correcting that node corrects one copy and leaves three that the table cannot reveal.

Searching the code as text finds all four.

> An hour spent this way is worth keeping.

### 2. Code That Computes Nothing

A node is a point where a quantity changes meaning.

Therefore, an error branch, a refusal condition, or a fallback has no node and appears nowhere in the table.

Such code can still change what is reported by leaving something out when it fires on real data.

Only **running the pipeline** shows this. That is Part A's job.

---

## 8 Granularity and Effort

### The Granularity Rule

A script of about **800 lines** usually yields **10–20 nodes**.

* Many more than that means the granularity is probably too fine. Merge neighbouring rows and let the regions grow.
* Many fewer means regions are absorbing decisions that should be visible as settings of their own.

The count is a **check on the enumeration, not a target**.

### Why Node Count Is Lower Than Operation Count

The node count is always lower than an operation count would suggest because the rules in Section 4.1 collapse repetition.

For example:

> One serial-correlation selection routine used by six scripts is one node and six call sites, not thirty separate rows.

### Effort

Enumeration is largely mechanical and runs at a few hours for a pipeline of twenty-odd scripts.

Checking, done forwards by someone who knows the data, runs at roughly **one to two minutes per node** once the table is in front of you. An afternoon therefore covers a table of a hundred rows across a few sittings.

Checking backwards takes several times that, which is the argument for:

> **Building in one direction and checking in the other.**

**Our own figures:** 14,000 lines across 22 scripts, *N* nodes. See note.

### Partial Application

A defensible partial application is to:

1. Enumerate everything.
2. Check exhaustively only the chains that reach reported quantities.
3. Mark the remainder as **enumerated but unchecked**.

That is still a stronger statement than any sample because it makes the unchecked portion explicit rather than invisible.
