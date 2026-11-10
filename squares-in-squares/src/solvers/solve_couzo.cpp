// Book solver: Francisco Couzo, certified best-bound witnesses (n=1..10000 index)
//   Source:  https://github.com/franciscouzo/square-packing (certificates/best, manifest.json)
//   License: MIT (repo LICENSE). Underlying packings carry the per-n credits in certificates/best/manifest.json (jlevy/squares, Ellsworth, Schadt, Muglich, PingYou Ltd, ...).
//   Data:    books/couzo.csv, fetched by tools/fetch_books.py on 2026-10-25. Exact rational certificates (x, y, t=tan(theta/2)) converted to float64 by tools/fetch_books.py; centres dilated by 1+1e-11 about the bounding-box middle to turn exact contacts into real gaps. Only the 720 counts with their own witness (n<=1000) are present; other n return no answer.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/couzo.csv", n); }
