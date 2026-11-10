// Book solver: Chris Deotte, "Squares in Squares 22.315 Data" - submission.csv (his 22.3152-level file, 2026-10-24)
//   Source:  https://www.kaggle.com/datasets/cdeotte/squares-in-squares-22-315-data   (see also his notebook "FTMC + fill_rect")
//   License: CC0-1.0 (dataset page).
//   Data:    books/deotte_22315.csv, fetched with the Kaggle CLI on 2026-10-26 (books/raw/deotte_22315/submission.csv).
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/deotte_22315.csv", n); }
