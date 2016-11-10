def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
def defaultCMakeOptions =
    '-DCMAKE_BUILD_TYPE=Release ' +
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System ' +
    '-DBUILD_SHARED_LIBS=ON'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

def image = docker.image('ogs6/gcc-gui:latest')
image.pull()
image.inside(defaultDockerArgs) {
    stage('Configure (Linux-Docker-Dynamic)') {
        configure.linux 'build', "${defaultCMakeOptions}"
    }

    stage('CLI (Linux-Docker-Dynamic)') {
        build.linux this, 'build'
    }

    stage('Test (Linux-Docker-Dynamic)') {
        build.linux this, 'build', 'tests ctest'
    }
}

stage('Post (Linux-Docker)') {
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules'
    post.cleanup()
}
