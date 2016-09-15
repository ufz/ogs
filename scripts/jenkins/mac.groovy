configure = load 'scripts/jenkins/lib/configure.groovy'
build     = load 'scripts/jenkins/lib/build.groovy'
post      = load 'scripts/jenkins/lib/post.groovy'
helper    = load 'scripts/jenkins/lib/helper.groovy'

defaultCMakeOptions = '-DCMAKE_BUILD_TYPE=Release -DOGS_CPU_ARCHITECTURE=core2 -DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System -DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'

node('mac && conan') {
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
        '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF',
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
