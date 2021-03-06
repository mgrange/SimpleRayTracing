solution "gKit2light"
	configurations { "debug", "release" }
	platforms { "x64", "x32" }

	includedirs { ".", "src/gKit" }

	configuration "debug"
		targetdir "bin/debug"
		defines { "DEBUG" }
		flags { "Symbols" }

	configuration "release"
		targetdir "bin/release"
--~ 		defines { "NDEBUG" }
--~ 		defines { "GK_RELEASE" }
		flags { "Optimize" }

	configuration "linux"
		buildoptions { "-mtune=native" }
		buildoptions { "-std=c++11" }
		buildoptions { "-W -Wall -Wextra -fopenmp  -Wsign-compare -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable", "-pipe" }
		buildoptions { "-flto"}
		linkoptions { "-flto -fopenmp"}
		links { "GLEW", "SDL2", "SDL2_image", "GL" }

	configuration { "linux", "debug" }
		buildoptions { "-g"}
		linkoptions { "-g"}

	configuration { "windows" }
		defines { "WIN32", "NVWIDGETS_EXPORTS", "_USE_MATH_DEFINES", "_CRT_SECURE_NO_WARNINGS" }
		defines { "NOMINMAX" } -- allow std::min() and std::max() in vc++ :(((

	configuration { "windows", "codeblocks", "x32" }
		includedirs { "extern/mingw/include" }
		libdirs { "extern/mingw/lib" }
		links { "mingw32", "SDL2main", "SDL2", "SDL2_image", "opengl32", "glew32" }

	configuration { "windows", "vs2013", "x64" }
		includedirs { "extern/visual2013/include" }
		libdirs { "extern/visual2013/lib" }
		links { "opengl32", "glew32", "SDL2", "SDL2main", "SDL2_image" }

	configuration { "windows", "vs2015", "x64" }
		includedirs { "extern/visual2015/include" }
		libdirs { "extern/visual2015/lib" }
		links { "opengl32", "glew32", "SDL2", "SDL2main", "SDL2_image" }

	configuration "macosx"
		local frameworks= "-F /Library/Frameworks/"
		defines { "GK_MACOS" }
		buildoptions { frameworks }
		linkoptions { frameworks .. " -framework OpenGL -framework SDL2 -framework SDL2_image" }


 -- description des fichiers communs
local gkit_files = { "src/gKit/*.cpp", "src/gKit/*.h" }

 -- description des RayTracing
local RayTracing = {
	"RayTracing"
}

for i, name in ipairs(RayTracing) do
	project(name)
		language "C++"
		kind "ConsoleApp"
		targetdir "bin"
		files ( gkit_files )
		files { "RayTracing/" .. name..'.cpp' }
		files { "RayTracing/" .. name..'Function.cpp' }
		files { "RayTracing/" .. name..'Struct.cpp' }
		files { "RayTracing/ReadingScene.cpp" }
end
