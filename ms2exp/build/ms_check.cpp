// Check one Orbit `ms query` output against MS truth at every position, streaming.
//
// Usage: ms_check OUTPUT TRUTH.u32 OFFSETS.u64
//
// OUTPUT has one line per read: its name q<i>, a tab, and its values separated by
// single spaces (Orbit QUERY_OUTPUT.md). TRUTH.u32 holds the little-endian uint32 MS
// of every read, concatenated in read order; OFFSETS.u64 holds the little-endian
// uint64 start of each read in TRUTH and the total, so read i spans
// [offsets[i], offsets[i + 1]). Every read must appear exactly once with exactly its
// length in values, each equal to truth. This is the check of real-hprc2
// benchmark.check for Orbit outputs, written to stream outputs with reads of hundreds
// of megabases (a whole chromosome) in constant memory.
//
// Prints one JSON object: {"reads":N,"positions":P,"failed":0} and exits 0 when
// everything matches; otherwise prints {"failed":1,"error":"..."} and exits 1.
#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& message) {
    std::printf("{\"failed\":1,\"error\":\"%s\"}\n", message.c_str());
    std::exit(1);
}

std::vector<uint64_t> read_offsets(const char* path) {
    FILE* f = std::fopen(path, "rb");
    if (!f) fail(std::string("cannot open offsets: ") + std::strerror(errno));
    std::vector<uint64_t> out;
    uint64_t v;
    while (std::fread(&v, sizeof v, 1, f) == 1) out.push_back(v);
    std::fclose(f);
    if (out.size() < 2 || out[0] != 0) fail("offsets must start at 0 and hold at least one read");
    for (size_t i = 1; i < out.size(); ++i)
        if (out[i] < out[i - 1]) fail("offsets decrease");
    return out;
}

// Buffered byte reader over a FILE*.
class Reader {
  public:
    explicit Reader(FILE* f) : f_(f), buf_(1 << 24) {}
    // Next byte, or -1 at end of file.
    int get() {
        if (pos_ == len_) {
            len_ = std::fread(buf_.data(), 1, buf_.size(), f_);
            pos_ = 0;
            if (len_ == 0) return -1;
        }
        return static_cast<unsigned char>(buf_[pos_++]);
    }

  private:
    FILE* f_;
    std::vector<char> buf_;
    size_t pos_ = 0, len_ = 0;
};

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: ms_check OUTPUT TRUTH.u32 OFFSETS.u64\n");
        return 2;
    }
    const std::vector<uint64_t> offsets = read_offsets(argv[3]);
    const size_t reads = offsets.size() - 1;
    FILE* truth = std::fopen(argv[2], "rb");
    if (!truth) fail(std::string("cannot open truth: ") + std::strerror(errno));
    if (std::fseek(truth, 0, SEEK_END) != 0) fail("cannot seek truth");
    if (static_cast<uint64_t>(std::ftell(truth)) != offsets.back() * 4) fail("truth size differs from the offsets");
    FILE* out = std::fopen(argv[1], "rb");
    if (!out) fail(std::string("cannot open output: ") + std::strerror(errno));
    Reader in(out);
    std::vector<char> seen(reads, 0);
    std::vector<uint32_t> expect(1 << 20);
    size_t lines = 0;
    for (;;) {
        int c = in.get();
        if (c == -1) break;
        // Name: q<i>, then a tab.
        if (c != 'q') fail("line " + std::to_string(lines) + ": name does not start with q");
        uint64_t i = 0;
        int digits = 0;
        while ((c = in.get()) >= '0' && c <= '9') {
            if (digits && i == 0) fail("line " + std::to_string(lines) + ": name has a leading zero");
            i = i * 10 + (c - '0');
            if (++digits > 19) fail("line " + std::to_string(lines) + ": name too long");
        }
        if (!digits || c != '\t') fail("line " + std::to_string(lines) + ": name is not q<i> followed by a tab");
        if (i >= reads) fail("unexpected read q" + std::to_string(i));
        if (seen[i]) fail("duplicate read q" + std::to_string(i));
        seen[i] = 1;
        const std::string name = "q" + std::to_string(i);
        const uint64_t length = offsets[i + 1] - offsets[i];
        if (std::fseek(truth, static_cast<long>(offsets[i] * 4), SEEK_SET) != 0) fail("cannot seek truth");
        uint64_t j = 0;          // values read so far
        size_t have = 0, at = 0;  // truth values buffered, next to compare
        bool more = true;
        while (more) {
            uint64_t value = 0;
            int vd = 0;
            while ((c = in.get()) >= '0' && c <= '9') {
                value = value * 10 + (c - '0');
                if (++vd > 10) fail(name + ": value too long");
            }
            if (!vd) fail(name + ": empty value at position " + std::to_string(j));
            if (c == '\n') more = false;
            else if (c != ' ') fail(name + ": unexpected character after position " + std::to_string(j));
            if (j == length) fail("length " + name + ": more values than its " + std::to_string(length) + " bases");
            if (at == have) {
                const uint64_t want = std::min<uint64_t>(expect.size(), length - j);
                if (std::fread(expect.data(), 4, want, truth) != want) fail("truth truncated");
                have = want;
                at = 0;
            }
            if (value != expect[at])
                fail("MS mismatch " + name + " at position " + std::to_string(j) + ": got " + std::to_string(value) +
                     ", truth " + std::to_string(expect[at]));
            ++at;
            ++j;
        }
        if (j != length) fail("length " + name + ": " + std::to_string(j) + " values for " + std::to_string(length) + " bases");
        ++lines;
    }
    std::fclose(out);
    std::fclose(truth);
    for (size_t i = 0; i < reads; ++i)
        if (!seen[i]) fail("missing reads: first missing q" + std::to_string(i));
    std::printf("{\"reads\":%zu,\"positions\":%llu,\"failed\":0}\n", reads,
                static_cast<unsigned long long>(offsets.back()));
    return 0;
}
