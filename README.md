# hs-diet

Harbour seal diet analysis project for the Skagerrak and nearby Swedish west coast areas.

This repository combines marine ecology, diet composition data, and spatial/spatiotemporal modelling to understand what harbour seals eat, where they eat it, and how diet varies with sex, size, time, and space.

## Project Focus

Main questions in this repository include:

- how harbour seal diet changes over time,
- how diet differs in space,
- how diet differs between prey families and prey taxa,
- whether seal sex, body length, and age relate to diet composition,
- whether seals are predicted to consume more cod, clupeids, or other prey in some areas, seasons, or life stages.

The analyses are primarily built around Quarto notebooks in `R/analyse-data/` and use `sdmTMB` for spatial and spatiotemporal models.

## Repository Structure

- `Data/`: source data, shapefiles, and imported datasets
- `R/analyse-data/`: main Quarto notebooks and analysis scripts
- `Literature/`: papers and supporting references
- `Presentation/`: presentation material
- `agent.md`: project-specific working guide for future edits
- `fixes.md`: running log of bugs, fixes, structural improvements, and modelling lessons

## Main Analysis Style

The preferred working style in this project is:

- marine-biological in interpretation,
- statistically cautious,
- explicit about assumptions,
- readable and stepwise,
- focused on clear model progression rather than dense compact code.

As a structural example, `03b-diet-models-family.qmd` is the clearest notebook pattern in the project and should be used as the main reference for organising new analyses.

`03c` contains useful analyses, but its structure became harder to follow. When extending the project, prefer the staged workflow style used in `03b`.

## Current Notebook Themes

Examples of current analysis directions include:

- family-level diet models,
- zero-one-inflated beta models for diet proportions,
- hunted-only ontogenetic models,
- hybrid prey-group models where Gadidae are partly split to scientific-name labels and the rest stay at family level,
- predictions on both relative diet composition and biomass-per-sample scales.

## Modelling Notes

The project frequently uses expanded sample-by-prey tables.

Important consequence:

- expanded modelling tables are correct for fitting models,
- but should not be used directly for seal-count summaries or histograms without `distinct(sample_id, ...)`.

For model interpretation, distinguish between:

- numerical convergence,
- residual fit,
- ecological plausibility.

Some models may converge well while still showing residual misfit. In those cases, use them for broad qualitative patterns rather than very precise effect-size claims.

## Logging Fixes And Lessons

This project should keep a running memory of what failed and what was improved.

Use:

- `fixes.md`

to record:

- bugs,
- failed model structures,
- formula simplifications,
- prediction fixes,
- plotting mistakes,
- notebook structure improvements,
- repeated pitfalls to avoid.

This is preferred over a formal `CHANGELOG.md`, because the main need here is not release tracking but analysis memory.

## Working Guidance

See:

- `agent.md`

for the project-specific working guide, including:

- expected role and style,
- preferred notebook structure,
- how and when to log fixes,
- expectations for future analysis edits.

## Practical Advice For Future Edits

Before editing a notebook:

1. read the relevant notebook section fully,
2. check `fixes.md` for known pitfalls,
3. prefer small, explicit changes,
4. preserve working code when adding new analyses,
5. if a new model or chunk fails, log it in `fixes.md`.

## Current State

The repository currently contains:

- family-level diet modelling workflows,
- ontogenetic diet analyses for hunted seals,
- hybrid prey-group models for cod, whiting, Gadidae, and non-gadid prey groups,
- prediction outputs for both composition and biomass per sample.

The project is actively exploratory, so model structures and figures may evolve as support diagnostics and biological questions become clearer.
