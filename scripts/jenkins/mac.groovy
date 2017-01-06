def defaultCMakeOptions =
    '-DCMAKE_BUILD_TYPE=Release ' +
    '-DOGS_CPU_ARCHITECTURE=core2 ' +
    '-DOGS_LIB_BOOST=System' +
    '-DOGS_LIB_VTK=System ' +
    '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ' +
    '-DCMAKE_OSX_DEPLOYMENT_TARGET="10.11" '

def guiCMakeOptions =
    '-DOGS_BUILD_GUI=ON ' +
    '-DOGS_BUILD_UTILS=ON ' +
    '-DOGS_BUILD_TESTS=OFF ' +
    '-DOGS_BUILD_METIS=ON '

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

stage('Configure (Mac)') {
    configure.linux(
        cmakeOptions: defaultCMakeOptions,
        generator: 'Ninja',
        script: this,
        useConan: true
    )
}

stage('CLI (Mac)') {
    build.linux(cmd: 'ninja', script: this)
}

stage('Test (Mac)') {
    build.linux(cmd: 'ninja', script: this, target: 'tests ctest')
}

stage('Data Explorer (Mac)') {
        configure.linux(
            cmakeOptions: defaultCMakeOptions + guiCMakeOptions,
            generator: 'Ninja',
            keepDir: true,
            script: this
        )
    build.linux(cmd: 'ninja', script: this)
}

stage('Archive (Mac)') {
    archiveArtifacts 'build/*.tar.gz,build/*.dmg,build/conaninfo.txt'
}

stage('Post (Mac)') {
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/msvc-log-parser.rules'
    post.cleanup()
}
