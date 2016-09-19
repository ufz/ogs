configure = load 'scripts/jenkins/lib/configure.groovy'
build     = load 'scripts/jenkins/lib/build.groovy'

node('win && conan') {
    def defaultCMakeOptions =
        '-DCMAKE_BUILD_TYPE=Release ' +
        '-DOGS_32_BIT=ON ' +
        '-DOGS_LIB_BOOST=System ' +
        '-DOGS_LIB_VTK=System ' +
        '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'

    ws {
        dir('ogs') { unstash 'source' }

        stage 'Data Explorer 32-bit (Win)'
        withEnv(helper.getEnv('x32')) {
            configure.win 'build-32', "${defaultCMakeOptions} " +
                '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF',
                'Ninja', '-u -s build_type=Release -s compiler="Visual ' +
                'Studio" -s compiler.version=12 -s arch=x86'
            build.win 'build-32'
        }
        archive 'build-32/*.zip'
    }
}
