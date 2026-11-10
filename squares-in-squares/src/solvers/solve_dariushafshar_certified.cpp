// Book solver: Dariush Afshar, Squares in Squares: Certified Packings 1-1000 - index_packings_dilated.csv
//   Source:  https://www.kaggle.com/datasets/dariushafshar/squares-in-squares-certified-packings
//   License: See dataset page (not verified here).
//   Data:    books/dariushafshar_certified.csv, fetched by tools/fetch_books.py on 2026-10-25. Dilated certified packings, one per n.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/dariushafshar_certified.csv", n); }
