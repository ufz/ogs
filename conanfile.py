## Appveyor cache version: 1
from conans import ConanFile, CMake

class OpenGeoSysConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"

    options = {"gui": [True, False]}
    default_options = \
        "gui=False", \
        "Boost:header_only=True", \
        "Qt:xmlpatterns=True"

    def requirements(self):
        self.requires("Boost/[>=1.56.0]@lasote/stable")
        self.requires("Eigen3/3.2.9@bilke/stable")
        self.requires("VTK/[>=7.1]@bilke/stable")

        if self.options.gui:
            self.requires("Shapelib/1.3.0@bilke/stable")
            self.requires("libgeotiff/1.4.2@bilke/stable")
            self.requires("Qt/5.6.2@bilke/testing")

    def imports(self):
        self.copy(pattern="*.dll", dst="bin", src="bin")
        self.copy(pattern="*.dylib*", dst="bin", src="lib")
        self.copy(pattern="*.so*", dst="bin", src="lib")
        self.copy(pattern="*.framework*", dst="bin", src="lib")
        self.copy(pattern="*.dll", dst="bin/platforms", src="plugins/platforms")
        self.copy(pattern="*.dylib*", dst="bin/platforms", src="plugins/platforms")

    def build(self):
        cmake = CMake(self.settings)
        self.run('cmake "%s" %s' % (self.conanfile_directory, cmake.command_line))
        self.run('cmake --build . %s' % cmake.build_config)
