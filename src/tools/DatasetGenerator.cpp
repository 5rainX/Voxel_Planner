#include "VoxelPlannerAPI.h"
#include "Module2_Morphology/VoxelIO.h"
#include "data/DatasetWriter.h"
#include "environment/SdfVolume.h"
#include "optimizer/SqpProjector.h"
#include "optimizer/Trajectory.h"

#include <algorithm>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string map_path;
    std::filesystem::path output_directory;
    int requested_samples = 0;
    int max_attempts = 0;
    float busbar_width = 0.0F;
    float busbar_thickness = 0.0F;
    int max_paths = 1;
    unsigned int seed = 5489U;
    double safe_margin = 1.0;
};

void printUsage(const char* executable) {
    std::cerr
        << "Usage: " << executable
        << " <map.txt_or_vtk> <output_dir> <sample_count>"
        << " <max_attempts> <width> <thickness>"
        << " [max_paths] [seed] [safe_margin]\n";
}

int parseInt(const char* value, const char* name) {
    try {
        std::size_t consumed = 0U;
        const int parsed = std::stoi(value, &consumed);
        if (value[consumed] != '\0') {
            throw std::invalid_argument("trailing characters");
        }
        return parsed;
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("Invalid integer: ") + name);
    }
}

float parseFloat(const char* value, const char* name) {
    try {
        std::size_t consumed = 0U;
        const float parsed = std::stof(value, &consumed);
        if (value[consumed] != '\0') {
            throw std::invalid_argument("trailing characters");
        }
        return parsed;
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("Invalid float: ") + name);
    }
}

double parseDouble(const char* value, const char* name) {
    try {
        std::size_t consumed = 0U;
        const double parsed = std::stod(value, &consumed);
        if (value[consumed] != '\0') {
            throw std::invalid_argument("trailing characters");
        }
        return parsed;
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("Invalid double: ") + name);
    }
}

Options parseOptions(int argc, char* argv[]) {
    if (argc < 7 || argc > 10) {
        throw std::invalid_argument("Wrong argument count.");
    }
    Options options;
    options.map_path = argv[1];
    options.output_directory = argv[2];
    options.requested_samples = parseInt(argv[3], "sample_count");
    options.max_attempts = parseInt(argv[4], "max_attempts");
    options.busbar_width = parseFloat(argv[5], "width");
    options.busbar_thickness = parseFloat(argv[6], "thickness");
    if (argc >= 8) {
        options.max_paths = parseInt(argv[7], "max_paths");
    }
    if (argc >= 9) {
        options.seed = static_cast<unsigned int>(parseInt(argv[8], "seed"));
    }
    if (argc >= 10) {
        options.safe_margin = parseDouble(argv[9], "safe_margin");
    }

    if (options.map_path.empty() ||
        options.requested_samples <= 0 ||
        options.max_attempts <= 0 ||
        options.busbar_width <= 0.0F ||
        options.busbar_thickness <= 0.0F ||
        options.max_paths <= 0 ||
        options.safe_margin < 0.0) {
        throw std::invalid_argument("Arguments must be positive.");
    }
    return options;
}

std::string sampleId(std::size_t index) {
    std::ostringstream stream;
    stream << "sample_" << std::setw(6) << std::setfill('0') << index;
    return stream.str();
}

void hashBytes(
    std::uint64_t& hash,
    const unsigned char* bytes,
    std::size_t count) noexcept {
    constexpr std::uint64_t kPrime = 1099511628211ULL;
    for (std::size_t index = 0U; index < count; ++index) {
        hash ^= static_cast<std::uint64_t>(bytes[index]);
        hash *= kPrime;
    }
}

template <typename T>
void hashValue(std::uint64_t& hash, const T& value) noexcept {
    hashBytes(
        hash,
        reinterpret_cast<const unsigned char*>(&value),
        sizeof(T));
}

void hashString(std::uint64_t& hash, const std::string& value) noexcept {
    hashBytes(
        hash,
        reinterpret_cast<const unsigned char*>(value.data()),
        value.size());
}

std::string hexHash(std::uint64_t hash) {
    std::ostringstream stream;
    stream << std::hex << std::setw(16) << std::setfill('0') << hash;
    return stream.str();
}

std::string makeMapHash(
    const Options& options,
    const VoxelGrid& grid) {
    std::uint64_t hash = 1469598103934665603ULL;
    hashString(
        hash,
        std::filesystem::absolute(options.map_path).generic_string());
    const std::uint32_t width = grid.width();
    const std::uint32_t height = grid.height();
    const std::uint32_t depth = grid.depth();
    hashValue(hash, width);
    hashValue(hash, height);
    hashValue(hash, depth);
    if (std::filesystem::exists(options.map_path)) {
        const std::uintmax_t size =
            std::filesystem::file_size(options.map_path);
        hashValue(hash, size);
    }
    return hexHash(hash);
}

std::string sdfRelativePath(const std::string& mapHash) {
    return std::string("sdf/") + mapHash + ".bin";
}

std::vector<voxel_planner::Point3D> collectFreeVoxels(
    const VoxelGrid& grid,
    const environment::GlobalSdf& sdf,
    double safeMargin) {
    std::vector<voxel_planner::Point3D> freeVoxels;
    freeVoxels.reserve(grid.voxelCount());
    for (std::uint32_t z = 0U; z < grid.depth(); ++z) {
        for (std::uint32_t y = 0U; y < grid.height(); ++y) {
            for (std::uint32_t x = 0U; x < grid.width(); ++x) {
                const std::size_t index = grid.index(x, y, z);
                if (grid.isRawObstacle(index)) {
                    continue;
                }
                const environment::Vec3 pos{
                    static_cast<double>(x),
                    static_cast<double>(y),
                    static_cast<double>(z)};
                if (sdf.getDistance(pos) >= safeMargin) {
                    freeVoxels.push_back({
                        static_cast<int>(x),
                        static_cast<int>(y),
                        static_cast<int>(z)});
                }
            }
        }
    }
    return freeVoxels;
}

double squaredDistance(
    const voxel_planner::Point3D& lhs,
    const voxel_planner::Point3D& rhs) {
    const double dx = static_cast<double>(lhs.x - rhs.x);
    const double dy = static_cast<double>(lhs.y - rhs.y);
    const double dz = static_cast<double>(lhs.z - rhs.z);
    return dx * dx + dy * dy + dz * dz;
}

std::vector<environment::Vec3> continuousPath(
    const optimizer::PoseTrajectory& trajectory) {
    std::vector<environment::Vec3> points;
    points.reserve(trajectory.size());
    for (const optimizer::TrajectoryNode& node : trajectory.nodes()) {
        points.push_back(node.position);
    }
    return points;
}

voxel_planner::EndpointPose endpointPoseFromNode(
    const optimizer::TrajectoryNode& node) {
    return {
        {node.normal.x, node.normal.y, node.normal.z},
        {node.tangent.x, node.tangent.y, node.tangent.z}};
}

} // namespace

int main(int argc, char* argv[]) {
    Options options;
    try {
        options = parseOptions(argc, argv);
    } catch (const std::exception& error) {
        printUsage(argv[0]);
        std::cerr << "Argument error: " << error.what() << "\n";
        return 2;
    }

    try {
        const VoxelGrid rawGrid =
            module2_morphology::VoxelIO::loadVoxelMap(options.map_path);
        const environment::GlobalSdf sdf(rawGrid);
        const std::string mapHash = makeMapHash(options, rawGrid);
        const std::string sdfBinRelativePath = sdfRelativePath(mapHash);
        const std::vector<voxel_planner::Point3D> freeVoxels =
            collectFreeVoxels(rawGrid, sdf, options.safe_margin);
        if (freeVoxels.size() < 2U) {
            throw std::runtime_error(
                "Not enough safe free voxels to sample endpoints.");
        }

        const voxel_planner::ProcessedMap map =
            voxel_planner::loadMap(
                options.map_path,
                options.busbar_width,
                options.busbar_thickness);
        const data::DatasetWriter writer(options.output_directory);
        std::mt19937 rng(options.seed);
        std::uniform_int_distribution<std::size_t> pick(
            0U,
            freeVoxels.size() - 1U);

        std::size_t accepted = 0U;
        std::size_t rejected = 0U;
        int attempts = 0;
        while (accepted < static_cast<std::size_t>(options.requested_samples) &&
               attempts < options.max_attempts) {
            ++attempts;
            try {
                const voxel_planner::Point3D start = freeVoxels[pick(rng)];
                const voxel_planner::Point3D goal = freeVoxels[pick(rng)];
                if (start == goal || squaredDistance(start, goal) < 25.0) {
                    ++rejected;
                    continue;
                }

                const auto plan = voxel_planner::findPaths(
                    map,
                    start,
                    goal,
                    options.max_paths);
                if (plan.first != voxel_planner::PlanStatus::OK ||
                    plan.second.empty() ||
                    plan.second.front().path.empty()) {
                    ++rejected;
                    continue;
                }

                const optimizer::ConfigurationSpaceRefiner refiner;
                const optimizer::PoseTrajectory initial =
                    refiner.fromVoxelPath(plan.second.front().path);
                optimizer::SqpProjectorOptions projectorOptions;
                projectorOptions.safe_margin = options.safe_margin;
                projectorOptions.busbar_width =
                    static_cast<double>(options.busbar_width);
                projectorOptions.busbar_thickness =
                    static_cast<double>(options.busbar_thickness);
                projectorOptions.max_curvature = 1.0 / std::min(
                    1.5 * static_cast<double>(options.busbar_thickness),
                    0.557143 * static_cast<double>(options.busbar_width));
                projectorOptions.max_iterations = 80;
                projectorOptions.initial_step = 0.08;
                const optimizer::SqpProjector projector(projectorOptions);
                const optimizer::SqpProjectorResult projected =
                    projector.project(initial, sdf);
                if (!projected.success ||
                    !projected.constraints_satisfied ||
                    projected.trajectory.size() < 2U) {
                    ++rejected;
                    continue;
                }

                const std::string id = sampleId(accepted);
                data::TrainingSample sample;
                sample.sample_id = id;
                sample.sdf_map_hash = mapHash;
                sample.sdf_bin_relative_path = sdfBinRelativePath;
                sample.sdf_origin = sdf.origin();
                sample.sdf_voxel_size = sdf.voxelSize();
                sample.sdf_width = sdf.dimensions()[0];
                sample.sdf_height = sdf.dimensions()[1];
                sample.sdf_depth = sdf.dimensions()[2];
                sample.start = start;
                sample.goal = goal;
                sample.start_pose =
                    endpointPoseFromNode(projected.trajectory.nodes().front());
                sample.goal_pose =
                    endpointPoseFromNode(projected.trajectory.nodes().back());
                sample.busbar_width = options.busbar_width;
                sample.busbar_thickness = options.busbar_thickness;
                sample.optimized_path = continuousPath(projected.trajectory);
                sample.voxel_path = refiner.toVoxelPath(projected.trajectory);
                const optimizer::PoseTrajectory discretizedTrajectory =
                    refiner.fromVoxelPath(sample.voxel_path, false);
                if (!projector.satisfiesHardConstraints(
                        discretizedTrajectory,
                        sdf)) {
                    ++rejected;
                    continue;
                }

                writer.writeSample(sample, sdf);
                ++accepted;
                std::cout << "[ACCEPT] " << id
                          << " attempt=" << attempts
                          << " path_nodes=" << sample.optimized_path.size()
                          << "\n";
            } catch (const std::exception& error) {
                ++rejected;
                std::cerr << "[REJECT] attempt=" << attempts
                          << " reason=\"" << error.what() << "\"\n";
                continue;
            } catch (...) {
                ++rejected;
                std::cerr << "[REJECT] attempt=" << attempts
                          << " reason=\"unknown exception\"\n";
                continue;
            }
        }

        std::cout << "Generated " << accepted << "/"
                  << options.requested_samples
                  << " samples after " << attempts
                  << " attempts; rejected=" << rejected << ".\n";
        return accepted == static_cast<std::size_t>(options.requested_samples)
            ? 0
            : 1;
    } catch (const std::exception& error) {
        std::cerr << "Dataset generation failed: " << error.what() << "\n";
        return 1;
    }
}
