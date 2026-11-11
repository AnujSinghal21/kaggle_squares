# Squares in Squares

A solver and analysis project for the Kaggle competition "Squares in Squares."

The task is to pack `n` unit squares into the smallest possible axis-aligned square while minimizing the score:

`L / sqrt(n) - 1`

where `L` is the side length of the bounding square. Lower is better.

This repository contains:
- C++ solver code for generating valid packings
- validation tools for checking submission legality
- multiple solver variants and public-book baselines
- Python utilities for analyzing scores and visualizing packings

## Repository layout

- `src/` — core solver and validation code
- `src/solvers/` — implementation variants and reference approaches
- `analysis/` — Python utilities for score plots and SVG visuals
- `tools/` — helpers for downloading public packings / books
- `data/` — sample Kaggle data
- `submissions/` — generated output submissions
- `Makefile` — build and execution entry points
- `books.mk` — list of book-based solver variants

## Getting started

From the repository root:

```bash
cd squares-in-squares
