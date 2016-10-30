node('docker') {
    def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
    def defaultCMakeOptions =
        '-DOGS_LIB_BOOST=System ' +
        '-DOGS_LIB_VTK=System ' +
        '-DOGS_BUILD_UTILS=ON ' +
        '-DOGS_NO_EXTERNAL_LIBS=ON'

    stage 'Checkout (Linux-Docker-No-Deps)'
    dir('ogs') { checkout scm }

    docker.image('ogs6/gcc-base:latest').inside(defaultDockerArgs) {
        stage 'Configure (Linux-Docker-No Deps)'
        configure.linux 'build', "${defaultCMakeOptions}"

        stage 'CLI (Linux-Docker-No-Deps)'
        build.linux 'build'

        stage 'Test (Linux-Docker-No-Deps)'
        build.linux 'build', 'tests ctest'
    }

    stage 'Post (Linux-Docker-No-Deps)'
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules'
    post.cleanup()
}
