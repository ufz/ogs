
def defaultCMakeOptions =
    '-DCMAKE_BUILD_TYPE=Release ' +
    '-DMSVC_RUNTIME=static ' +
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System ' +
    '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'

def guiCMakeOptions =
    '-DOGS_BUILD_GUI=ON ' +
    '-DOGS_BUILD_UTILS=ON ' +
    '-DOGS_BUILD_TESTS=OFF ' +
    '-DOGS_BUILD_SWMM=ON ' +
    '-DOGS_BUILD_METIS=ON'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

withEnv(helper.getEnv(this)) {
    stage('Configure (Win)') {
        configure.win 'build', "${defaultCMakeOptions}", 'Ninja',
            '-u -s build_type=Release -s compiler="Visual Studio" ' +
            '-s compiler.version=12 -s arch=x86_64'
    }

    stage('CLI (Win)') {
        build.win this, 'build'
    }

    stage('Test (Win)') {
        build.win this, 'build', 'tests'
        build.win this, 'build', 'ctest'
    }

    stage('Data Explorer (Win)') {
        configure.win 'build', "${defaultCMakeOptions} ${guiCMakeOptions}",
            'Ninja', '-u -s build_type=Release -s compiler="Visual Studio" ' +
            '-s compiler.version=12 -s arch=x86_64', true
        build.win this, 'build'
    }
}

stage('Archive (Win)') {
    archiveArtifacts 'build/*.zip'
}

stage('Post (Win)') {
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/msvc-log-parser.rules'
    post.cleanup()
}
