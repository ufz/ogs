#include "../partition.h"

namespace MeshLib::IO
{
bool isFileManager()
{
    return true;
}

std::pair<std::size_t, std::size_t> getPartitionInfo(std::size_t const size)
{
    return {0, size};
}
}  // namespace MeshLib::IO