defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
defaultCMakeOptions = '-DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System'

node('docker') {
    def build = new ogs.build()

    stage('Checkout') { dir('ogs') { checkout scm } }

    stage('Build') {
        docker.image('ogs6/gcc-base:latest').inside(defaultDockerArgs) {
            build this, 'build', '-DOGS_COVERAGE=ON',
                'testrunner_coverage_cobertura ctest_coverage_cobertura'
        }
    }

    archiveArtifacts 'build/*.xml'
    // TODO: Report is published in a free-style child job
}
