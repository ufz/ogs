def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache ' +
    '-v /home/jenkins/conan-data:/home/jenkins/.conan/data'

def defaultCMakeOptions =
    '-DCMAKE_BUILD_TYPE=Release ' +
    '-DOGS_USE_CONAN=ON ' +
    '-DOGS_CONAN_BUILD=never ' +
    '-DOGS_CPU_ARCHITECTURE=generic ' +
    '-DOGS_PACKAGE_DEPENDENCIES=ON '

def guiCMakeOptions =
    '-DOGS_BUILD_CLI=OFF ' +
    '-DOGS_BUILD_GUI=ON ' +
    '-DOGS_BUILD_UTILS=ON ' +
    '-DOGS_BUILD_TESTS=OFF ' +
    '-DOGS_BUILD_METIS=ON '

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

def image = docker.image('ogs6/gcc-conan')
image.pull()
image.inside(defaultDockerArgs) {
    stage('Configure (Linux-Docker)') {
        configure.linux(
            cmakeOptions: defaultCMakeOptions,
            script: this
        )
    }

    stage('CLI (Linux-Docker)') {
        build.linux(script: this)
    }

    stage('Test (Linux-Docker)') {
        build.linux(script: this, target: 'tests ctest')
    }

    stage('Data Explorer (Linux-Docker)') {
        configure.linux(
            cmakeOptions: defaultCMakeOptions + guiCMakeOptions,
            keepDir: true,
            script: this
        )
        build.linux(script: this)
    }
}

stage('Archive (Linux-Docker)') {
    archiveArtifacts 'build/*.tar.gz,build/conaninfo.txt'
}

stage('Post (Linux-Docker)') {
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml'
    post.cleanup()
}
