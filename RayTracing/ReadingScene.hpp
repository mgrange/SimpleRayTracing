#ifndef _READING_SCENE_H
#define _READING_SCENE_H

#include <string>

#include "RayTracingStruct.hpp"
#include "RayTracingFunction.hpp"

void readScene(const std::string filename, listscene & allelement, light_and_cam & lumiere);

#endif
