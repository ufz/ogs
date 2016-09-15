configure = load 'scripts/jenkins/lib/configure.groovy'
build     = load 'scripts/jenkins/lib/build.groovy'
post      = load 'scripts/jenkins/lib/post.groovy'
helper    = load 'scripts/jenkins/lib/helper.groovy'

defaultCMakeOptions = '-DCMAKE_BUILD_TYPE=Release -DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System ' +
    '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'

node('win && conan') {
    stage 'Checkout (Win)'
    dir('ogs') { checkout scm }

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
        configure.win 'build', "${defaultCMakeOptions} " +
            '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF',
            'Ninja', '-u -s build_type=Release -s compiler="Visual Studio" -s compiler.version=12' +
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
