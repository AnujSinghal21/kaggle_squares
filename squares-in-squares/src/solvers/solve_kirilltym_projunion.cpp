// Book solver: Kirill Tymchishin (kirilltym), Squares in Squares | CPU/CUDA/MPS | EDA + ProjUnion (notebook output)
//   Source:  https://www.kaggle.com/code/kirilltym/squares-in-squares-cpu-cuda-mps-eda-projunion
//   License: See notebook page and its output licenses/LICENSE-LEVY-SOURCES.txt.
//   Data:    books/kirilltym_projunion.csv, fetched by tools/fetch_books.py on 2026-10-25. Projected union of public witnesses plus CPU contact-LP refinement.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/kirilltym_projunion.csv", n); }
