// A minimal tar.gz writer, for archives Python's tarfile writes with
// tarfile.open(path, "w:gz") (chicExportData).
//
// Members are regular files in the POSIX ustar layout; a name longer than the
// ustar fields hold gets a pax extended header with a path record, as
// tarfile's default PAX_FORMAT writes it. The archive ends with two zero
// blocks and is padded to tarfile's record size of 20 blocks, then gzip
// compressed at level 9. Member metadata (mode, owner, modification time) is
// what the caller passes; nothing compares it, because tarfile stamps
// time.time().

#ifndef HICX_TAR_GZ_HPP
#define HICX_TAR_GZ_HPP

#include <cstdint>
#include <string>
#include <string_view>

namespace hicx {

class TarGzWriter {
  public:
    explicit TarGzWriter(const std::string& path);
    ~TarGzWriter();
    TarGzWriter(const TarGzWriter&) = delete;
    TarGzWriter& operator=(const TarGzWriter&) = delete;

    // A regular file member with the given content.
    void add(const std::string& name, std::string_view content, std::int64_t mtime,
             unsigned mode = 0644);
    // A regular file member with the content of a file on disk.
    void add_file(const std::string& name, const std::string& source, std::int64_t mtime,
                  unsigned mode = 0644);
    // Writes the end of archive blocks and closes the file; throws on error.
    void close();

  private:
    void write(const char* data, std::size_t length);
    void header(const std::string& name, std::uint64_t size, std::int64_t mtime, unsigned mode,
                char type);

    std::string path_;
    void* gz_ = nullptr;
    std::uint64_t written_ = 0;
};

}  // namespace hicx

#endif  // HICX_TAR_GZ_HPP
