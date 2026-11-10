// Book solver: alex chilton, Squares in Squares starter kit - best.csv
//   Source:  https://www.kaggle.com/datasets/alexchilton/squares-in-squares-starter-kit
//   License: See dataset page (not verified here). Builds on Couzo (MIT), Evan Daniel (MIT) and Levy (MIT code / CC BY 4.0 data).
//   Data:    books/alexchilton_starter_kit.csv, fetched by tools/fetch_books.py on 2026-10-25. Complete 1..1000 submission file.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/alexchilton_starter_kit.csv", n); }
