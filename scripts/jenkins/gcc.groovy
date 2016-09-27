node('docker') {
    def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
    def defaultCMakeOptions =
        '-DOGS_LIB_BOOST=System ' +
        '-DOGS_LIB_VTK=System'

    stage 'Checkout (Linux-Docker)'
    dir('ogs') { checkout scm }

    docker.image('ogs6/gcc-gui:latest').inside(defaultDockerArgs) {
        stage 'Configure (Linux-Docker)'
        configure.linux 'build', "${defaultCMakeOptions}"

        stage 'CLI (Linux-Docker)'
        build.linux 'build'

        stage 'Test (Linux-Docker)'
        build.linux 'build', 'tests ctest'

        stage 'Data Explorer (Linux-Docker)'
        configure.linux 'build', "${defaultCMakeOptions} " +
            '-DOGS_BUILD_CLI=OFF -DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON ' +
            '-DOGS_BUILD_TESTS=OFF -DOGS_BUILD_METIS=ON',
            'Unix Makefiles', null, true
        build.linux 'build'
    }

    if (helper.isRelease()) {
        stage 'Release (Linux-Docker)'
        archive 'build/*.tar.gz'
    }

    stage 'Post (Linux-Docker)'
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules'
    post.cleanup()
}
