from conans import ConanFile, CMake

class OpenGeoSysConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"

    requires = \
        "Boost/[>=1.56.0]@lasote/stable", \
        "Shapelib/1.3.0@bilke/stable", \
        "VTK/[>=6.3,<7.1]@bilke/stable", \
        "Eigen3/3.2.9@bilke/stable", \
        "libgeotiff/1.4.2@bilke/stable", \
        "Qt/5.6.2@bilke/testing"

    generators = "cmake"

    default_options = \
        "Boost:header_only=True", \
        "Qt:xmlpatterns=True"

    def imports(self):
        self.copy(pattern="*.dll", dst="bin", src="bin")
        self.copy(pattern="*.dylib*", dst="bin", src="lib")
        # self.copy(pattern="*.framework*", dst="bin", src="lib")
        self.copy(pattern="*.dll", dst="bin/platforms", src="plugins/platforms")
        self.copy(pattern="*.dylib*", dst="bin/platforms", src="plugins/platforms")

    def build(self):
        cmake = CMake(self.settings)
        self.run('cmake "%s" %s' % (self.conanfile_directory, cmake.command_line))
        self.run('cmake --build . %s' % cmake.build_config)
