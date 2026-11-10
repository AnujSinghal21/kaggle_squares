// Book solver: prvsiyan, Squares in Squares Auto-updater live state - submission.csv
//   Source:  https://www.kaggle.com/datasets/prvsiyan/squares-auto-updater-live
//   License: See dataset page (not verified here). Per-size credits are in that dataset's credits.csv.
//   Data:    books/prvsiyan_auto_updater.csv, fetched by tools/fetch_books.py on 2026-10-25. Live-updated public portfolio; the fetched copy is a snapshot from 2026-10-25.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/prvsiyan_auto_updater.csv", n); }
