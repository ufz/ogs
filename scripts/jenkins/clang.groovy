node('docker') {
    def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
    def defaultCMakeOptions =
        '-DOGS_LIB_BOOST=System ' +
        '-DOGS_LIB_VTK=System ' +
        '-DOGS_ADDRESS_SANITIZER=ON ' +
        '-DOGS_UNDEFINED_BEHAVIOR_SANITIZER=ON'

    stage 'Checkout (Clang)'
    dir('ogs') { checkout scm }

    docker.image('ogs6/clang-base:latest').inside(defaultDockerArgs) {
        stage 'Configure (Clang)'
        configure.linux 'build', "${defaultCMakeOptions}"
        try {
            stage 'Unit tests (Clang)'
            build.linux 'build', 'tests', 'UBSAN_OPTIONS=print_stacktrace=1 make -j $(nproc)'
        }
        catch(err) {
            echo "Clang sanitizer for unit tests failed, marking build as unstable!"
            currentBuild.result = "UNSTABLE"
        }

        try {
            stage 'End-to-end tests (Clang)'
            build.linux 'build', 'ctest', 'UBSAN_OPTIONS=print_stacktrace=1 make -j $(nproc)'
        }
        catch(err) {
            echo "Clang sanitizer for end-to-end tests failed, marking build as unstable!"
            currentBuild.result = "UNSTABLE"
        }

    }

    stage 'Post (Clang)'
    post.publishTestReports('build/Testing/**/*.xml','build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules')
}
