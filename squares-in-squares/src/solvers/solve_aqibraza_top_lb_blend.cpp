// Book solver: Aqib Raza, Squares in Squares - Top LB Blending + Visualization (notebook output)
//   Source:  https://www.kaggle.com/code/aqibrazadev/squares-in-squares-top-lb-blending-visualization
//   License: Notebook output; Kaggle notebooks are Apache-2.0 unless stated. Inherits credits of the public notebooks it blends.
//   Data:    books/aqibraza_top_lb_blend.csv, fetched by tools/fetch_books.py on 2026-10-25. Blend of public submissions, ~22.2988 on the public leaderboard when fetched.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/aqibraza_top_lb_blend.csv", n); }
