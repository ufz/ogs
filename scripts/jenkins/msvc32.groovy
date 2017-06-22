def defaultCMakeOptions =
    '-DCMAKE_BUILD_TYPE=Release ' +
    '-DOGS_32_BIT=ON ' +
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System ' +
    '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ' +
    '-DOGS_BUILD_GUI=ON ' +
    '-DOGS_BUILD_UTILS=OFF ' +
    '-DOGS_BUILD_TESTS=OFF ' +
    '-DOGS_BUILD_CLI=OFF'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()


withEnv(helper.getEnv(this, 'x32')) {
    stage('Data Explorer 32-bit (Win)') {
        configure.win(
            arch: 'x86',
            cmakeOptions: defaultCMakeOptions,
            conanOptions: "-o gui=True",
            script: this
        )
        build.win(script: this)
    }
}

stage('Post 32-bit (Win)') {
    archiveArtifacts 'build/*.zip'
    post.cleanup('build')
}
