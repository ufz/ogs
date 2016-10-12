node('docker') {
    stage 'Checkout (Clang)'
    dir('ogs') {
        checkout scm

        configure = load 'scripts/jenkins/lib/configure.groovy'
        build     = load 'scripts/jenkins/lib/build.groovy'
        post      = load 'scripts/jenkins/lib/post.groovy'
    }

    def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
    def defaultCMakeOptions =
        '-DOGS_LIB_BOOST=System ' +
        '-DOGS_LIB_VTK=System ' +
        '-DOGS_ADDRESS_SANITIZER=ON ' +
        '-DOGS_UNDEFINED_BEHAVIOR_SANITIZER=ON ' +
        '-DOGS_BUILD_UTILS=ON'

    docker.image('ogs6/clang-base:latest').inside(defaultDockerArgs) {
        stage 'Configure (Clang)'
        configure.linux 'build', "${defaultCMakeOptions}"
        try {
            stage 'Unit tests (Clang)'
            build.linux 'build', 'tests', 'UBSAN_OPTIONS=print_stacktrace=1 make -j $(nproc)'
        }
        catch(err) { echo "Clang sanitizer for unit tests failed!" }

        try {
            stage 'End-to-end tests (Clang)'
            build.linux 'build', 'ctest', 'UBSAN_OPTIONS=print_stacktrace=1 make -j $(nproc)'
        }
        catch(err) { echo "Clang sanitizer for end-to-end tests failed!" }
    }

    stage 'Post (Clang)'
    post.publishTestReports('build/Testing/**/*.xml','build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules')
    post.cleanup()
}

properties([[
    $class: 'org.jenkinsci.plugins.workflow.job.properties.BuildDiscarderProperty',
    strategy: [$class: 'LogRotator',
    artifactDaysToKeepStr: '',
    artifactNumToKeepStr: '5',
    daysToKeepStr: '',
    numToKeepStr: '25']]])
