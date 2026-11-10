// Book solver: lightningv08, Squares in Squares Layout Atlas - alex_merge.csv
//   Source:  https://www.kaggle.com/datasets/lightningv08/squares-in-squares-layout-atlas
//   License: See dataset page (the author's CUDA search code is Apache-2.0).
//   Data:    books/lightningv08_atlas_merge.csv, fetched by tools/fetch_books.py on 2026-10-25. Merge of public packings with provenance in alex_merge_sources.csv.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/lightningv08_atlas_merge.csv", n); }
