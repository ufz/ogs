#!/usr/bin/env groovy

node('docker') {
    def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
    def defaultCMakeOptions =
        '-DOGS_LIB_BOOST=System ' +
        '-DOGS_LIB_VTK=System ' +
        '-DOGS_ADDRESS_SANITIZER=ON ' +
        '-DOGS_UNDEFINED_BEHAVIOR_SANITIZER=ON ' +
        '-DOGS_BUILD_UTILS=ON'

    def configure = new ogs.configure()
    def build = new ogs.build()
    def post = new ogs.post()
    def helper = new ogs.helper()

    stage('Checkout (Clang)') {
        dir('ogs') { checkout scm }
    }

    def image = docker.image('ogs6/clang-base:latest')
    image.pull()
    image.inside(defaultDockerArgs) {
        stage('Configure (Clang)') {
            configure.linux 'build', "${defaultCMakeOptions}"
        }
        try {
            stage('Unit tests (Clang)') {
                build.linux this, 'build', 'tests', 'UBSAN_OPTIONS=print_stacktrace=1 make -j $(nproc)'
            }
        }
        catch(err) { echo "Clang sanitizer for unit tests failed!" }

        try {
            stage('End-to-end tests (Clang)') {
                build.linux this, 'build', 'ctest', 'UBSAN_OPTIONS=print_stacktrace=1 make -j $(nproc)'
            }
        }
        catch(err) { echo "Clang sanitizer for end-to-end tests failed!" }
    }

    stage('Post (Clang)') {
        post.publishTestReports('build/Testing/**/*.xml','build/Tests/testrunner.xml',
            'ogs/scripts/jenkins/clang-log-parser.rules')
        post.cleanup()
    }
}
