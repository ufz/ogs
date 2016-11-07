def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
def defaultCMakeOptions =
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System ' +
    '-DOGS_USE_LIS=ON'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

node('docker') {
    def image = docker.image('ogs6/gcc-gui:latest')
    image.pull()
    image.inside(defaultDockerArgs) {
        stage('Configure') { configure.linux 'build', "${defaultCMakeOptions}" }
        stage('Build') { build.linux this, 'build' }
        stage('Test') { build.linux this, 'build', 'tests ctest-large' }
    }

    stage('Post') {
        post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
            'ogs/scripts/jenkins/clang-log-parser.rules'
        post.cleanup()
    }
}
