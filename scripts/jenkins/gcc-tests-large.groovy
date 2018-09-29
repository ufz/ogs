def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
def defaultCMakeOptions =
    '-DCMAKE_BUILD_TYPE=Release ' +
    '-DOGS_USE_LIS=ON'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

node('envinf11w') {
    checkout scm
    def image = docker.image('ogs6/gcc-gui:latest')
    image.pull()
    image.inside(defaultDockerArgs) {
        stage('Configure') { configure.linux(cmakeOptions: defaultCMakeOptions,
                                             script: this) }
        stage('Build') { build.linux(script: this) }
        stage('Test') { build.linux(cmd: 'make -j 1',
                                    script: this,
                                    target: 'tests ctest-large-serial') }
    }

    stage('Post') {
        post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml'
        post.cleanup()
    }
}
