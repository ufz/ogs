def defaultCMakeOptions =
    '-DCMAKE_BUILD_TYPE=Release ' +
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System ' +
    '-DBUILD_SHARED_LIBS=ON ' +
    '-DOGS_BUILD_UTILS=ON'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

stage('Configure (Linux-Docker-Dynamic)') {
    configure.linux(cmakeOptions: defaultCMakeOptions, env: 'envinf1/cli.sh', script: this)
}

stage('CLI (Linux-Docker-Dynamic)') {
    build.linux(env: 'envinf1/cli.sh', script: this)
}

stage('Test (Linux-Docker-Dynamic)') {
    build.linux(env: 'envinf1/cli.sh', script: this, target: 'tests ctest')
}

stage('Post (Linux-Docker-Dynamic)') {
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules'
    post.cleanup()
}
