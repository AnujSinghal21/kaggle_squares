// Book solver: Vladislav Lavrentev (lightningv08), Squares in Squares | CUDA Packing Search (notebook output)
//   Source:  https://www.kaggle.com/code/lightningv08/squares-in-squares-cuda-packing-search
//   License: Apache-2.0 (stated by the author); background credit to prvsiyan.
//   Data:    books/lightningv08_cuda_search.csv, fetched by tools/fetch_books.py on 2026-10-25. CUDA L-BFGS search plus fault-line recomposition.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/lightningv08_cuda_search.csv", n); }
