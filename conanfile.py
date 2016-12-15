from conans import ConanFile, CMake

class OpenGeoSysConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"

    requires = \
        "Boost/[>=1.56.0]@lasote/stable", \
        "Shapelib/1.3.0@bilke/stable", \
        "VTK/[>=6.3,<7.1]@bilke/stable", \
        "Eigen3/3.2.8@bilke/stable", \
        "libgeotiff/1.4.1@bilke/stable"

    generators = "cmake"

    default_options = \
        "Boost:shared=False", \
        "Boost:header_only=True"

    def build(self):
        cmake = CMake(self.settings)
        self.run('cmake "%s" %s' % (self.conanfile_directory, cmake.command_line))
        self.run('cmake --build . %s' % cmake.build_config)
