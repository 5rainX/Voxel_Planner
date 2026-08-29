#pragma once

#include "Module1_Main/Types.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <ios>
#include <limits>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

namespace module2_morphology {

/**
 * @brief Owns voxel file parsing and VTK serialization.
 *
 * Raw storage values remain private to VoxelGrid and this trusted I/O class.
 */
class VoxelIO {
public:
    static VoxelGrid loadVoxelMap(const std::string& inputPath) {
        const std::string extension = lowercaseExtension(inputPath);
        if (extension == ".vtk") {
            return loadAsciiVtk(inputPath);
        }
        return loadTextVoxelMap(inputPath);
    }

    static void exportDilatedMapToVtk(
        const VoxelGrid& grid,
        const std::string& outputPath) {
        std::string header;
        header.reserve(256);
        header += "# vtk DataFile Version 3.0\nDual-kernel voxel map\n"
                  "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS ";
        appendInteger(header, grid.width_);
        header.push_back(' ');
        appendInteger(header, grid.height_);
        header.push_back(' ');
        appendInteger(header, grid.depth_);
        header += "\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA ";
        appendInteger(header, grid.storage_.size());
        header += "\nSCALARS voxel unsigned_char 1\n"
                  "LOOKUP_TABLE default\n";

        std::string output(
            header.size() + grid.storage_.size() * 2U,
            '\0');
        std::memcpy(output.data(), header.data(), header.size());

        const std::size_t baseOffset = header.size();
        const std::int64_t voxelCount =
            static_cast<std::int64_t>(grid.storage_.size());
#pragma omp parallel for schedule(static)
        for (std::int64_t i = 0; i < voxelCount; ++i) {
            const std::size_t offset =
                baseOffset + static_cast<std::size_t>(i) * 2U;
            const std::uint8_t value =
                grid.storage_[static_cast<std::size_t>(i)];
            output[offset] = value >= 3U ? '1' : '0';
            output[offset + 1U] =
                ((i + 1) % 16 == 0 || i + 1 == voxelCount) ? '\n' : ' ';
        }

        std::ofstream out(outputPath, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            throw std::runtime_error(
                "Failed to open VTK output file: " + outputPath);
        }
        out.write(output.data(), static_cast<std::streamsize>(output.size()));
        if (!out) {
            throw std::runtime_error(
                "Failed while writing VTK output file: " + outputPath);
        }
    }

    static void exportVoxelClassificationToVtk(
        const VoxelGrid& grid,
        const std::string& outputPath) {
        std::string header;
        header.reserve(256);
        header += "# vtk DataFile Version 3.0\n"
                  "Pose-aware voxel classification\n"
                  "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS ";
        appendInteger(header, grid.width_);
        header.push_back(' ');
        appendInteger(header, grid.height_);
        header.push_back(' ');
        appendInteger(header, grid.depth_);
        header += "\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA ";
        appendInteger(header, grid.storage_.size());
        header += "\nSCALARS VoxelClass int 1\n"
                  "LOOKUP_TABLE default\n";

        std::vector<int> classifications(grid.storage_.size(), 0);
        for (std::uint32_t z = 0; z < grid.depth_; ++z) {
            for (std::uint32_t y = 0; y < grid.height_; ++y) {
                for (std::uint32_t x = 0; x < grid.width_; ++x) {
                    const std::size_t linearIndex = grid.index(x, y, z);
                    classifications[linearIndex] =
                        voxelClassToScalar(grid.getState(linearIndex));
                }
            }
        }

        std::ofstream out(outputPath, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            throw std::runtime_error(
                "Failed to open VTK output file: " + outputPath);
        }
        out << header;
        for (std::size_t i = 0; i < classifications.size(); ++i) {
            out << classifications[i];
            out << ((i + 1U) % 16U == 0U ||
                    i + 1U == classifications.size()
                        ? '\n'
                        : ' ');
        }
        if (!out) {
            throw std::runtime_error(
                "Failed while writing VTK output file: " + outputPath);
        }
    }

    static void exportVoxelClassificationAndPosesToVtk(
        const VoxelGrid& grid,
        const std::string& outputPath) {
        std::string header;
        header.reserve(320);
        header += "# vtk DataFile Version 3.0\n"
                  "Pose-aware voxel classification and pose counts\n"
                  "ASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS ";
        appendInteger(header, grid.width_);
        header.push_back(' ');
        appendInteger(header, grid.height_);
        header.push_back(' ');
        appendInteger(header, grid.depth_);
        header += "\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA ";
        appendInteger(header, grid.storage_.size());
        header += "\nSCALARS VoxelClass int 1\n"
                  "LOOKUP_TABLE default\n";

        std::vector<int> classifications(grid.storage_.size(), 0);
        std::vector<int> allowedPoseCounts(grid.storage_.size(), 0);
        for (std::uint32_t z = 0; z < grid.depth_; ++z) {
            for (std::uint32_t y = 0; y < grid.height_; ++y) {
                for (std::uint32_t x = 0; x < grid.width_; ++x) {
                    const std::size_t linearIndex = grid.index(x, y, z);
                    const VoxelState state = grid.getState(linearIndex);
                    classifications[linearIndex] =
                        voxelClassToScalar(state);
                    if (state == VoxelState::UNCONDITIONAL) {
                        allowedPoseCounts[linearIndex] =
                            static_cast<int>(grid.poseCount_);
                    } else if (state == VoxelState::POSE_CONDITIONAL) {
                        allowedPoseCounts[linearIndex] =
                            static_cast<int>(
                                grid.allowedPoseCount(linearIndex));
                    }
                }
            }
        }

        std::ofstream out(outputPath, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            throw std::runtime_error(
                "Failed to open VTK output file: " + outputPath);
        }
        out << header;
        writeScalarPayload(out, classifications);
        out << "\nSCALARS AllowedPoseCount int 1\n"
               "LOOKUP_TABLE default\n";
        writeScalarPayload(out, allowedPoseCounts);
        if (!out) {
            throw std::runtime_error(
                "Failed while writing VTK output file: " + outputPath);
        }
    }

private:
    class MappedTextFile {
    public:
        explicit MappedTextFile(const std::string& inputPath) {
#ifdef _WIN32
            file_ = CreateFileA(
                inputPath.c_str(),
                GENERIC_READ,
                FILE_SHARE_READ,
                nullptr,
                OPEN_EXISTING,
                FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN,
                nullptr);
            if (file_ == INVALID_HANDLE_VALUE) {
                throw std::runtime_error(
                    "Failed to open text voxel map: " + inputPath);
            }

            LARGE_INTEGER fileSize{};
            if (!GetFileSizeEx(file_, &fileSize) || fileSize.QuadPart <= 0 ||
                static_cast<unsigned long long>(fileSize.QuadPart) >
                    static_cast<unsigned long long>(
                        std::numeric_limits<std::size_t>::max())) {
                closeHandles();
                throw std::runtime_error(
                    "Failed to determine text voxel map size.");
            }
            size_ = static_cast<std::size_t>(fileSize.QuadPart);

            mapping_ = CreateFileMappingA(
                file_, nullptr, PAGE_READONLY, 0, 0, nullptr);
            if (mapping_ == nullptr) {
                closeHandles();
                throw std::runtime_error(
                    "Failed to create read-only file mapping.");
            }
            view_ = static_cast<const char*>(MapViewOfFile(
                mapping_, FILE_MAP_READ, 0, 0, 0));
            if (view_ == nullptr) {
                closeHandles();
                throw std::runtime_error(
                    "Failed to map text voxel map into memory.");
            }
#else
            std::ifstream in(inputPath, std::ios::binary | std::ios::ate);
            if (!in.is_open()) {
                throw std::runtime_error(
                    "Failed to open text voxel map: " + inputPath);
            }
            const std::streampos endPosition = in.tellg();
            if (endPosition <= 0) {
                throw std::runtime_error(
                    "Failed to determine text voxel map size.");
            }
            fallback_.resize(static_cast<std::size_t>(endPosition));
            in.seekg(0, std::ios::beg);
            in.read(fallback_.data(),
                    static_cast<std::streamsize>(fallback_.size()));
            if (!in) {
                throw std::runtime_error(
                    "Failed to read text voxel map: " + inputPath);
            }
            view_ = fallback_.data();
            size_ = fallback_.size();
#endif
        }

        ~MappedTextFile() {
#ifdef _WIN32
            closeHandles();
#endif
        }

        MappedTextFile(const MappedTextFile&) = delete;
        MappedTextFile& operator=(const MappedTextFile&) = delete;

        const char* data() const noexcept { return view_; }
        std::size_t size() const noexcept { return size_; }

    private:
#ifdef _WIN32
        void closeHandles() noexcept {
            if (view_ != nullptr) {
                UnmapViewOfFile(view_);
                view_ = nullptr;
            }
            if (mapping_ != nullptr) {
                CloseHandle(mapping_);
                mapping_ = nullptr;
            }
            if (file_ != INVALID_HANDLE_VALUE) {
                CloseHandle(file_);
                file_ = INVALID_HANDLE_VALUE;
            }
        }

        HANDLE file_ = INVALID_HANDLE_VALUE;
        HANDLE mapping_ = nullptr;
#endif
        const char* view_ = nullptr;
        std::size_t size_ = 0;
#ifndef _WIN32
        std::string fallback_;
#endif
    };

    static std::string lowercaseExtension(const std::string& path) {
        const std::size_t dot = path.find_last_of('.');
        if (dot == std::string::npos) {
            return {};
        }
        std::string extension = path.substr(dot);
        std::transform(
            extension.begin(),
            extension.end(),
            extension.begin(),
            [](unsigned char character) {
                return static_cast<char>(std::tolower(character));
            });
        return extension;
    }

    static void skipWhitespace(const char*& cursor, const char* end) {
        while (cursor < end &&
               static_cast<unsigned char>(*cursor) <= ' ') {
            ++cursor;
        }
    }

    static bool parseInteger(
        const char*& cursor,
        const char* end,
        int& value) {
        skipWhitespace(cursor, end);
        if (cursor >= end) {
            return false;
        }

        int sign = 1;
        if (*cursor == '-' || *cursor == '+') {
            sign = *cursor == '-' ? -1 : 1;
            ++cursor;
        }
        if (cursor >= end || *cursor < '0' || *cursor > '9') {
            return false;
        }

        int parsed = 0;
        while (cursor < end && *cursor >= '0' && *cursor <= '9') {
            parsed = parsed * 10 + (*cursor - '0');
            ++cursor;
        }
        value = sign * parsed;
        return true;
    }

    static std::string_view parseToken(
        const char*& cursor,
        const char* end) {
        skipWhitespace(cursor, end);
        const char* begin = cursor;
        while (cursor < end &&
               static_cast<unsigned char>(*cursor) > ' ') {
            ++cursor;
        }
        return std::string_view(
            begin,
            static_cast<std::size_t>(cursor - begin));
    }

    static VoxelGrid loadTextVoxelMap(const std::string& inputPath) {
        std::ifstream input(inputPath);
        if (!input.is_open()) {
            throw std::runtime_error(
                "Failed to open text voxel map: " + inputPath);
        }

        std::string headerLine;
        if (!std::getline(input, headerLine)) {
            throw std::runtime_error("Text voxel map is empty.");
        }
        if (headerLine.size() >= 3U &&
            static_cast<unsigned char>(headerLine[0]) == 0xEFU &&
            static_cast<unsigned char>(headerLine[1]) == 0xBBU &&
            static_cast<unsigned char>(headerLine[2]) == 0xBFU) {
            headerLine.erase(0U, 3U);
        }

        // Both "map: x y z" and "header x y z" are accepted. The first
        // token is a descriptive prefix; only the three dimensions matter.
        std::stringstream header(headerLine);
        std::string prefix;
        int width = 0;
        int height = 0;
        int depth = 0;
        if (!(header >> prefix >> width >> height >> depth)) {
            throw std::runtime_error(
                "Expected map dimensions after the header prefix.");
        }
        if (width <= 0 || height <= 0 || depth <= 0) {
            throw std::runtime_error("Text voxel map has invalid dimensions.");
        }

        VoxelGrid grid(
            static_cast<std::uint32_t>(width),
            static_cast<std::uint32_t>(height),
            static_cast<std::uint32_t>(depth));

        std::string line;
        while (std::getline(input, line)) {
            std::stringstream record(line);
            int x = 0;
            int y = 0;
            int z = 0;
            if (!(record >> x >> y >> z)) {
                continue;
            }
            if (!grid.isValid(x, y, z)) {
                continue;
            }

            // Sparse text records are occupancy declarations. Any trailing
            // value is metadata and must not affect occupancy.
            grid.storage_[grid.index(
                static_cast<std::uint32_t>(x),
                static_cast<std::uint32_t>(y),
                static_cast<std::uint32_t>(z))] =
                VoxelGrid::kRawObstacleValue;
        }
        grid.rebuildOccupancyBlocks();
        return grid;
    }

    static VoxelGrid loadAsciiVtk(const std::string& inputPath) {
        std::ifstream in(inputPath);
        if (!in.is_open()) {
            throw std::runtime_error(
                "Failed to open VTK voxel map: " + inputPath);
        }

        std::string line;
        std::uint32_t width = 0;
        std::uint32_t height = 0;
        std::uint32_t depth = 0;
        bool ascii = false;
        bool datasetSupported = false;
        bool scalarPayloadFound = false;
        while (std::getline(in, line)) {
            std::istringstream stream(line);
            std::string keyword;
            if (!(stream >> keyword)) {
                continue;
            }
            if (keyword == "ASCII") {
                ascii = true;
            } else if (keyword == "DATASET") {
                std::string type;
                stream >> type;
                datasetSupported =
                    type == "STRUCTURED_POINTS" ||
                    type == "ImageData" ||
                    type == "IMAGE_DATA";
            } else if (keyword == "DIMENSIONS") {
                stream >> width >> height >> depth;
            } else if (keyword == "LOOKUP_TABLE") {
                scalarPayloadFound = true;
                break;
            }
        }

        if (!ascii || !datasetSupported || width == 0 || height == 0 ||
            depth == 0 || !scalarPayloadFound) {
            throw std::runtime_error(
                "Unsupported VTK format. Expected ASCII structured volume.");
        }

        VoxelGrid grid(width, height, depth);
        for (std::size_t i = 0; i < grid.storage_.size(); ++i) {
            int scalar = 0;
            if (!(in >> scalar)) {
                throw std::runtime_error("VTK scalar payload ended early.");
            }
            grid.storage_[i] = scalar != 0 ? 255U : 0U;
        }
        grid.rebuildOccupancyBlocks();
        return grid;
    }

    template <typename Integer>
    static void appendInteger(std::string& output, Integer value) {
        char buffer[32];
        const auto converted = std::to_chars(
            buffer,
            buffer + sizeof(buffer),
            value);
        if (converted.ec != std::errc()) {
            throw std::runtime_error("Failed to format VTK integer.");
        }
        output.append(buffer, converted.ptr);
    }

    static int voxelClassToScalar(VoxelState state) noexcept {
        switch (state) {
        case VoxelState::BLOCKED:
            return 0;
        case VoxelState::UNCONDITIONAL:
            return 1;
        case VoxelState::POSE_CONDITIONAL:
            return 2;
        }
        return 0;
    }

    static void writeScalarPayload(
        std::ofstream& out,
        const std::vector<int>& values) {
        for (std::size_t i = 0; i < values.size(); ++i) {
            out << values[i];
            out << ((i + 1U) % 16U == 0U ||
                    i + 1U == values.size()
                        ? '\n'
                        : ' ');
        }
    }
};

} // namespace module2_morphology
