// Book solver: alex chilton, Squares in Squares starter kit - quench.csv: puzzles improved by quenching fault-DP gadgets
//   Source:  https://www.kaggle.com/datasets/alexchilton/squares-in-squares-starter-kit
//   License: dataset page says "other" (code CC0; packings are the cited authors'); local competition use only, not verified for redistribution.
//   Data:    books/starterkit_quench.csv, fetched with the Kaggle CLI on 2026-10-26 (books/raw/starter_kit_full/quench.csv).
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/starterkit_quench.csv", n); }
