#include "../../include/utils.h"

#include <sys/stat.h>   // For mkdir on POSIX
#include <sys/types.h>  // For mode_t
#include <errno.h>
#include <string>
#include <vector>
#include <cstring>      // For strerror
#include <iostream>     // For debug printing

#ifdef _WIN32
  #include <direct.h>   // For _mkdir
  #define MKDIR(name, mode) _mkdir(name)
#else
  #include <unistd.h>   // For access, etc.
  #define MKDIR(name, mode) mkdir(name, mode)
#endif

/**
 * @brief Create all directories in the given path (like std::filesystem::create_directories).
 * Created 2025 by ChatGPT to emulate a capability found in C++17.
 * @param dirPath The directory path, using forward slashes or backslashes on Windows.
 * @return true if successful or directory already existed, false on error.
 */
bool create_directories(const std::string &dirPath)
{
    // An empty string or a single slash doesn't need to do anything
    if (dirPath.empty()) {
        return true;
    }

    // Choose a directory separator. On Windows, both '/' and '\\' can work.
    // We'll just look for either slash so we can unify splitting logic.
    // In practice, you might normalize your path first.
    char sep = '/';
#ifdef _WIN32
    // Could also check for backslash if you want to handle paths with '\'
    // But forward slashes generally work, so we’ll keep it simple.
#endif

    // We'll split on the separator to handle each folder in the hierarchy
    std::vector<std::string> parts;
    {
        // Replace backslashes with forward slashes to unify
        std::string normalized;
        normalized.reserve(dirPath.size());
        for (char c : dirPath) {
            if (c == '\\')
                normalized.push_back('/');
            else
                normalized.push_back(c);
        }

        size_t start = 0;
        for (size_t i = 0; i < normalized.size(); i++) {
            if (normalized[i] == sep) {
                // Avoid empty tokens (could happen if multiple slashes)
                if (i > start) {
                    parts.push_back(normalized.substr(start, i - start));
                }
                start = i + 1;
            }
        }
        // Last token
        if (start < normalized.size()) {
            parts.push_back(normalized.substr(start));
        }
    }

    // For absolute paths on Unix, the first part might be empty.
    // For Windows drive letters (e.g., "C:"), you might need special handling.
    // We'll keep it simple here.

    // Accumulate subdirectories one at a time
    std::string partial;
    for (size_t i = 0; i < parts.size(); i++) {
        if (!partial.empty()) {
            partial += sep;
        }
        partial += parts[i];

        // Attempt to create this subfolder
        if (MKDIR(partial.c_str(), 0755) != 0) {
            // If already exists, ignore the error
            // On Windows: EEXIST is 17, same as POSIX
            if (errno != EEXIST) {
                // Another error?
                // e.g., EACCES means permission denied
                // EROFS means read-only filesystem, etc.
                std::cerr << "Error creating directory `" << partial 
                          << "`: " << strerror(errno) << std::endl;
                return false;
            }
        }
    }

    return true;
}