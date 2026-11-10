// Book solver: Sung Heo (kmmmvs), Squares Dashboard Public Batch 20261025 (notebook output)
//   Source:  https://www.kaggle.com/code/kmmmvs/squares-dashboard-public-batch-20261025
//   License: Apache-2.0 plus Levy source notices (see the companion dataset kmmmvs/squares-in-squares-best-packings: LICENSE-APACHE-2.0.txt, LICENSE-LEVY-SOURCES.txt).
//   Data:    books/kmmmvs_dashboard.csv, fetched by tools/fetch_books.py on 2026-10-25. Aggregates public best-validated packings with per-contributor credit.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/kmmmvs_dashboard.csv", n); }
