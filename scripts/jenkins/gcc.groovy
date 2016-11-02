def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
def defaultCMakeOptions =
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

docker.image('ogs6/gcc-gui:latest').inside(defaultDockerArgs) {
    stage('Configure (Linux-Docker)') {
        configure.linux 'build', "${defaultCMakeOptions}"
    }

    stage('CLI (Linux-Docker)') {
        build.linux this, 'build'
    }

    stage('Test (Linux-Docker)') {
        build.linux this, 'build', 'tests ctest'
    }

    stage('Data Explorer (Linux-Docker)') {
        configure.linux 'build', "${defaultCMakeOptions} " +
            '-DOGS_BUILD_CLI=OFF -DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON ' +
            '-DOGS_BUILD_TESTS=OFF -DOGS_BUILD_METIS=ON',
            'Unix Makefiles', null, true
        build.linux this, 'build'
    }
}

if (helper.isRelease(this)) {
    stage('Release (Linux-Docker)') {
        archiveArtifacts 'build/*.tar.gz'
    }
}

stage('Post (Linux-Docker)') {
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules'
    post.cleanup()
}
