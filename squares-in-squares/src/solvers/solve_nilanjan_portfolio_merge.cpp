// Book solver: Nilanjan Saha, S in S: Certified Portfolio Merge + Proj Quench (notebook output)
//   Source:  https://www.kaggle.com/code/nilanjansaha123/s-in-s-certified-portfolio-merge-proj-quench
//   License: Notebook output; see notebook page. Inherits credits of merged public packings.
//   Data:    books/nilanjan_portfolio_merge.csv, fetched by tools/fetch_books.py on 2026-10-25. Contributes whole-puzzle starts for n=545 and n=869 per kirilltym's notes.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/nilanjan_portfolio_merge.csv", n); }
