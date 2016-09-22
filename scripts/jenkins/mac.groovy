node('mac && conan') {
    def defaultCMakeOptions =
        '-DCMAKE_BUILD_TYPE=Release ' +
        '-DOGS_CPU_ARCHITECTURE=core2 ' +
        '-DOGS_LIB_BOOST=System' +
        '-DOGS_LIB_VTK=System ' +
        '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ' +
        '-DCMAKE_OSX_DEPLOYMENT_TARGET="10.11"'

    stage 'Checkout (Mac)'
    dir('ogs') { checkout scm }

    stage 'Configure (Mac)'
    configure.linux 'build', "${defaultCMakeOptions}", 'Ninja', ''

    stage 'CLI (Mac)'
    build.linux 'build', null, 'ninja'

    stage 'Test (Mac)'
    build.linux 'build', 'tests ctest', 'ninja'

    stage 'Data Explorer (Mac)'
    configure.linux 'build', "${defaultCMakeOptions} " +
        '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF ' +
        '-DOGS_BUILD_METIS=ON',
        'Ninja', '', true
    build.linux 'build', null, 'ninja'

    if (helper.isRelease()) {
        stage 'Release (Mac)'
        archive 'build/*.tar.gz,build/*.dmg'
    }

    stage 'Post (Mac)'
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/msvc-log-parser.rules'
}
