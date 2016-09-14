configure = load 'scripts/jenkins/lib/configure.groovy'
build     = load 'scripts/jenkins/lib/build.groovy'
post      = load 'scripts/jenkins/lib/post.groovy'
helper    = load 'scripts/jenkins/lib/helper.groovy'

defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
defaultCMakeOptions = '-DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System'

node('docker') {
    stage 'Checkout (Linux-Docker)'
    dir('ogs') { checkout scm }

    docker.image('ogs6/gcc-base:latest').inside(defaultDockerArgs) {
        stage 'Configure (Linux-Docker)'
        configure.linux 'build', ''

        stage 'CLI (Linux-Docker)'
        build.linux 'build', ''

        stage 'Test (Linux-Docker)'
        build.linux 'build', 'tests ctest'

        if (env.BRANCH_NAME == 'master' || env.BRANCH_NAME.contains('release')) {
            stage 'Release (Linux-Docker)'
            build.linux 'build', 'package'
        }
    }

    if (helper.isRelease())
        archive 'build/*.tar.gz'

    stage 'Post (Linux-Docker)'
    helper.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules'
}
