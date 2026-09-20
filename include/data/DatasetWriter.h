#pragma once

#include "VoxelPlannerAPI.h"
#include "environment/SdfVolume.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace data {

struct TrainingSample {
    std::string sample_id;
    std::string sdf_map_hash;
    std::string sdf_bin_relative_path;
    environment::Vec3 sdf_origin{};
    environment::Vec3 sdf_voxel_size{1.0, 1.0, 1.0};
    std::size_t sdf_width = 0U;
    std::size_t sdf_height = 0U;
    std::size_t sdf_depth = 0U;

    voxel_planner::Point3D start{};
    voxel_planner::Point3D goal{};
    voxel_planner::EndpointPose start_pose{};
    voxel_planner::EndpointPose goal_pose{};
    float busbar_width = 0.0F;
    float busbar_thickness = 0.0F;

    std::vector<environment::Vec3> optimized_path;
    std::vector<voxel_planner::Point3D> voxel_path;
};

class DatasetWriter {
public:
    explicit DatasetWriter(std::filesystem::path output_directory)
        : output_directory_(std::move(output_directory)) {
        std::filesystem::create_directories(output_directory_);
    }

    void writeSample(
        const TrainingSample& sample,
        const environment::GlobalSdf& sdf) const {
        if (sample.sample_id.empty()) {
            throw std::invalid_argument("TrainingSample sample_id is empty.");
        }
        if (sample.optimized_path.empty() || sample.voxel_path.empty()) {
            throw std::invalid_argument(
                "TrainingSample path payload must not be empty.");
        }

        TrainingSample normalized = sample;
        if (normalized.sdf_bin_relative_path.empty()) {
            normalized.sdf_bin_relative_path =
                makeSdfRelativePath(normalized.sdf_map_hash);
        }
        const std::filesystem::path binPath =
            output_directory_ / normalized.sdf_bin_relative_path;
        const std::filesystem::path jsonPath =
            output_directory_ / (normalized.sample_id + ".json");
        std::filesystem::create_directories(binPath.parent_path());
        std::filesystem::create_directories(jsonPath.parent_path());

        writeSdfBinaryIfMissing(binPath, sdf.data());
        writeJson(jsonPath, normalized);
    }

private:
    static std::string makeSdfRelativePath(const std::string& mapHash) {
        if (mapHash.empty()) {
            throw std::invalid_argument(
                "TrainingSample sdf_map_hash is empty.");
        }
        return std::string("sdf/") + jsonEscape(mapHash) + ".bin";
    }

    static std::string nativeEndianness() noexcept {
        const std::uint16_t value = 1U;
        const auto* bytes =
            reinterpret_cast<const unsigned char*>(&value);
        return bytes[0] == 1U ? "little" : "big";
    }

    static std::string jsonEscape(const std::string& value) {
        std::ostringstream escaped;
        for (const char ch : value) {
            switch (ch) {
            case '\\':
                escaped << "\\\\";
                break;
            case '"':
                escaped << "\\\"";
                break;
            case '\n':
                escaped << "\\n";
                break;
            case '\r':
                escaped << "\\r";
                break;
            case '\t':
                escaped << "\\t";
                break;
            default:
                escaped << ch;
                break;
            }
        }
        return escaped.str();
    }

    static std::uintmax_t expectedBinarySize(
        const std::vector<float>& values) {
        if (values.size() >
            static_cast<std::size_t>(
                std::numeric_limits<std::streamsize>::max()) /
                sizeof(float)) {
            throw std::overflow_error(
                "SDF binary payload is too large to write.");
        }
        return static_cast<std::uintmax_t>(
            values.size() * sizeof(float));
    }

    static std::size_t checkedElementCount(const TrainingSample& sample) {
        const auto checkedMultiply = [](
                                     std::size_t lhs,
                                     std::size_t rhs) {
            if (lhs != 0U &&
                rhs > std::numeric_limits<std::size_t>::max() / lhs) {
                throw std::overflow_error(
                    "TrainingSample SDF dimensions overflow.");
            }
            return lhs * rhs;
        };
        return checkedMultiply(
            checkedMultiply(sample.sdf_width, sample.sdf_height),
            sample.sdf_depth);
    }

    static void writeSdfBinaryIfMissing(
        const std::filesystem::path& path,
        const std::vector<float>& values) {
        const std::uintmax_t expectedSize = expectedBinarySize(values);
        if (std::filesystem::exists(path)) {
            const std::uintmax_t existingSize =
                std::filesystem::file_size(path);
            if (existingSize != expectedSize) {
                throw std::runtime_error(
                    "Existing SDF binary has unexpected size: " +
                    path.string());
            }
            return;
        }

        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            throw std::runtime_error(
                "Failed to open SDF binary output: " + path.string());
        }
        out.write(
            reinterpret_cast<const char*>(values.data()),
            static_cast<std::streamsize>(expectedSize));
        if (!out) {
            throw std::runtime_error(
                "Failed while writing SDF binary output: " + path.string());
        }
    }

    static void writeVec3Json(
        std::ofstream& out,
        const environment::Vec3& value) {
        out << '[' << value.x << ", " << value.y << ", " << value.z << ']';
    }

    static void writePointJson(
        std::ofstream& out,
        const voxel_planner::Point3D& value) {
        out << '[' << value.x << ", " << value.y << ", " << value.z << ']';
    }

    static void writePoseJson(
        std::ofstream& out,
        const voxel_planner::EndpointPose& pose) {
        out << "{\n      \"normal\": ";
        writeVec3Json(out, {pose.normal.x, pose.normal.y, pose.normal.z});
        out << ",\n      \"tangent\": ";
        writeVec3Json(out, {pose.tangent.x, pose.tangent.y, pose.tangent.z});
        out << "\n    }";
    }

    static void writeJson(
        const std::filesystem::path& path,
        const TrainingSample& sample) {
        std::ofstream out(path, std::ios::trunc);
        if (!out.is_open()) {
            throw std::runtime_error(
                "Failed to open JSON output: " + path.string());
        }

        out << std::fixed << std::setprecision(8);
        out << "{\n";
        out << "  \"sample_id\": \"" << jsonEscape(sample.sample_id)
            << "\",\n";
        out << "  \"sdf\": {\n";
        out << "    \"map_hash\": \"" << jsonEscape(sample.sdf_map_hash)
            << "\",\n";
        out << "    \"bin\": \""
            << jsonEscape(sample.sdf_bin_relative_path) << "\",\n";
        out << "    \"dtype\": \"float32\",\n";
        out << "    \"endianness\": \"" << nativeEndianness() << "\",\n";
        out << "    \"layout\": \"z_y_x_flat\",\n";
        out << "    \"origin\": ";
        writeVec3Json(out, sample.sdf_origin);
        out << ",\n    \"voxel_size\": ";
        writeVec3Json(out, sample.sdf_voxel_size);
        out << ",\n    \"W\": " << sample.sdf_width
            << ",\n    \"H\": " << sample.sdf_height
            << ",\n    \"D\": " << sample.sdf_depth
            << ",\n    \"element_count\": "
            << checkedElementCount(sample)
            << "\n";
        out << "  },\n";
        out << "  \"task\": {\n";
        out << "    \"busbar_width\": " << sample.busbar_width
            << ",\n    \"busbar_thickness\": "
            << sample.busbar_thickness << ",\n";
        out << "    \"start\": ";
        writePointJson(out, sample.start);
        out << ",\n    \"goal\": ";
        writePointJson(out, sample.goal);
        out << ",\n    \"start_pose\": ";
        writePoseJson(out, sample.start_pose);
        out << ",\n    \"goal_pose\": ";
        writePoseJson(out, sample.goal_pose);
        out << "\n  },\n";
        out << "  \"path\": {\n";
        out << "    \"voxel\": [\n";
        for (std::size_t i = 0U; i < sample.voxel_path.size(); ++i) {
            out << "      ";
            writePointJson(out, sample.voxel_path[i]);
            out << (i + 1U == sample.voxel_path.size() ? "\n" : ",\n");
        }
        out << "    ],\n";
        out << "    \"continuous\": [\n";
        for (std::size_t i = 0U; i < sample.optimized_path.size(); ++i) {
            out << "      ";
            writeVec3Json(out, sample.optimized_path[i]);
            out << (i + 1U == sample.optimized_path.size() ? "\n" : ",\n");
        }
        out << "    ]\n";
        out << "  }\n";
        out << "}\n";
        if (!out) {
            throw std::runtime_error(
                "Failed while writing JSON output: " + path.string());
        }
    }

    std::filesystem::path output_directory_;
};

} // namespace data
