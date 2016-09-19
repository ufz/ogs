node('win && conan') {
    def defaultCMakeOptions =
        '-DCMAKE_BUILD_TYPE=Release ' +
        '-DOGS_LIB_BOOST=System ' +
        '-DOGS_LIB_VTK=System ' +
        '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'

    def guiCMakeOptions =
        '-DOGS_BUILD_GUI=ON ' +
        '-DOGS_BUILD_UTILS=ON ' +
        '-DOGS_BUILD_TESTS=OFF ' +
        '-DOGS_BUILD_SWMM=ON'

    ws {
        dir('ogs') { unstash 'source' }

        withEnv(helper.getEnv()) {
            stage 'Configure (Win)'
            configure.win 'build', "${defaultCMakeOptions}", 'Ninja',
                '-u -s build_type=Release -s compiler="Visual Studio" -s compiler.version=12 -s ' +
                    'arch=x86_64'

            stage 'CLI (Win)'
            build.win 'build'

            stage 'Test (Win)'
            build.win 'build', 'tests'

            stage 'Data Explorer (Win)'
            configure.win 'build', "${defaultCMakeOptions} ${guiCMakeOptions}", 'Ninja',
                '-u -s build_type=Release -s compiler="Visual Studio" -s compiler.version=12' +
                ' -s arch=x86_64',
                true
            build.win 'build'
        }

        if (helper.isRelease()) {
            stage 'Release (Win)'
            archive 'build/*.zip'
        }

        stage 'Post (Win)'
        post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
            'ogs/scripts/jenkins/msvc-log-parser.rules'
    }
}
