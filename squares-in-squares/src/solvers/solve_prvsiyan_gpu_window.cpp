// Book solver: prvsiyan, A Square Deal - Frozen Public Inputs - gpu_window_v2.csv
//   Source:  https://www.kaggle.com/datasets/prvsiyan/a-square-deal-frozen-public-inputs
//   License: See dataset NOTICE.txt on the dataset page.
//   Data:    books/prvsiyan_gpu_window.csv, fetched by tools/fetch_books.py on 2026-10-25. GPU window-search results; possibly partial, missing n return no answer.
// Returns the book's packing for n, or an empty packing if the book has none. Validity is checked
// by main.cpp before anything is adopted, so a bad entry can never worsen a submission.
#include "../book.hpp"
#include "../solve.hpp"

Packing solve(int n, const Packing&) { return book_packing("books/prvsiyan_gpu_window.csv", n); }
