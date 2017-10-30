
def defaultCMakeOptions =
    '-DOGS_USE_CONAN=ON ' +
    '-DOGS_CONAN_BUILD=never ' +
    '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ' +
    '-DOGS_PACKAGE_DEPENDENCIES=ON '

def guiCMakeOptions =
    '-DOGS_BUILD_GUI=ON ' +
    '-DOGS_BUILD_UTILS=ON ' +
    '-DOGS_BUILD_TESTS=OFF ' +
    '-DOGS_BUILD_SWMM=ON ' +
    '-DOGS_BUILD_METIS=ON '

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

withEnv(helper.getEnv(this)) {
    stage('Configure (Win)') {
        configure.win(cmakeOptions: defaultCMakeOptions, script: this)
    }

    stage('CLI (Win)') {
        build.win(script: this)
    }

    stage('Test (Win)') {
        build.win(script: this, target: 'tests')
        build.win(script: this, target: 'ctest')
    }

    stage('Data Explorer (Win)') {
        configure.win(
            cmakeOptions: defaultCMakeOptions + guiCMakeOptions, keepDir: true,
            script: this
        )
        build.win(script: this)
    }
}

stage('Archive (Win)') {
    archiveArtifacts 'build/*.zip,build/conaninfo.txt'
}

stage('Post (Win)') {
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml'
    post.cleanup()
}
