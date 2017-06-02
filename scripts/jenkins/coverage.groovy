defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'

def configure = new ogs.configure()
def build = new ogs.build()

def image = docker.image('ogs6/gcc-base:latest')
image.pull()
image.inside(defaultDockerArgs) {
    stage('Configure (Coverage)') {
        configure.linux(cmakeOptions: '-DOGS_COVERAGE=ON', script: this)
    }

    stage('Build (Coverage)') {
        build.linux(
            script: this,
            target: 'testrunner_coverage_cobertura ctest_coverage_cobertura'
        )
    }
}

stage('Publish (Coverage)') {
    step([$class: 'CoberturaPublisher', coberturaReportFile: 'build/*.xml'])
}
