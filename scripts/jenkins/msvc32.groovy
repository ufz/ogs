def defaultCMakeOptions =
    '-DCMAKE_BUILD_TYPE=Release ' +
    '-DOGS_32_BIT=ON ' +
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System ' +
    '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()


withEnv(helper.getEnv(this, 'x32')) {
    stage('Data Explorer 32-bit (Win)') {
        configure.win 'build-32', "${defaultCMakeOptions} " +
            '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF ' +
            '-DOGS_BUILD_METIS=ON',
            'Ninja', '-u -s build_type=Release -s compiler="Visual ' +
            'Studio" -s compiler.version=12 -s arch=x86'
        build.win this, 'build-32'
    }
}

stage('Post 32-bit (Win)') {
    archiveArtifacts 'build-32/*.zip'
    post.cleanup('build-32')
}
