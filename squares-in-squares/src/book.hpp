// Shared helper for the "book" solvers: look a packing up in a normalised book CSV
// (books/<name>.csv, built by tools/fetch_books.py). The CSV is parsed once per process.
#pragma once
#include <string>

#include "common.hpp"
#include "validate.hpp"

// Packing of n squares from the book, or empty if the book has no complete packing for n.
inline Packing book_packing(const std::string& csv_path, int n) {
    static const Submission book = [&] {
        Submission s = load_submission(csv_path);
        if (!s.fatal.empty()) std::fprintf(stderr, "book %s: %s (run tools/fetch_books.py)\n", csv_path.c_str(), s.fatal.c_str());
        return s;
    }();
    if (!book.fatal.empty() || !book.present[n] || !book.broken[n].empty()) return {};
    return book.pk[n];
}
