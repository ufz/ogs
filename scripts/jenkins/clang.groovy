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
        catchError {
            configure.linux 'build', "${defaultCMakeOptions}"

            stage 'Unit tests (Clang)'
            build.linux 'build', 'tests', 'UBSAN_OPTIONS=print_stacktrace=1 make -j $(nproc)'

            stage 'End-to-end tests (Clang)'
            build.linux 'build', 'ctest', 'UBSAN_OPTIONS=print_stacktrace=1 make -j $(nproc)'
        }
    }

    stage 'Post (Clang)'
    post.publishTestReports('build/Testing/**/*.xml','build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules')
}
