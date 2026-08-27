#include "Module2_Morphology/VoxelIO.h"

#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        return 2;
    }

    const VoxelGrid grid =
        module2_morphology::VoxelIO::loadVoxelMap(argv[1]);
    if (grid.width() != 4U || grid.height() != 2U || grid.depth() != 1U ||
        grid.getState(0, 0, 0) != VoxelState::BLOCKED ||
        grid.getState(1, 0, 0) != VoxelState::BLOCKED ||
        grid.getState(2, 0, 0) != VoxelState::BLOCKED ||
        grid.getState(3, 0, 0) != VoxelState::UNCONDITIONAL) {
        std::cerr << "Missing-flag parser verification failed.\n";
        return 1;
    }

    std::cout << "[PASS] missing flag, flag=0, and flag=1 are obstacles.\n";
    return 0;
}
