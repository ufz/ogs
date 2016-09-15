configure = load 'scripts/jenkins/lib/configure.groovy'
build     = load 'scripts/jenkins/lib/build.groovy'
post      = load 'scripts/jenkins/lib/post.groovy'

defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'

node('docker') {
    stage 'Checkout'
    dir('ogs') { checkout scm }

    docker.image('ogs6/clang-base:latest').inside(defaultDockerArgs) {
        catchError {
            configure.linux 'build', "${defaultCMakeOptions} " +
            '-DOGS_ADDRESS_SANITIZER=ON -DOGS_UNDEFINED_BEHAVIOR_SANITIZER=ON',
            ''

            stage 'Unit tests'
            build.linux 'build', 'tests', 'UBSAN_OPTIONS=print_stacktrace=1 make'

            stage 'End-to-end tests'
            build.linux 'build', 'ctest', 'UBSAN_OPTIONS=print_stacktrace=1 make'
        }
    }
    post.publishTestReports('build/Testing/**/*.xml','build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules')
}
