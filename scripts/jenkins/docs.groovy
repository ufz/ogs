def defaultCMakeOptions =
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

docker.image('ogs6/gcc-base:16.04').inside() {
    stage('Configure (Docs)') {
        configure.linux 'build', "${defaultCMakeOptions}"
    }

    stage('Generate (Docs)') {
        build.linux 'build', 'doc'
    }
}

stage('Reports (Docs)') {
    publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: false, keepAll: false,
        reportDir: 'build/docs', reportFiles: 'index.html', reportName: 'Doxygen'])
    step([$class: 'WarningsPublisher', canResolveRelativePaths: false,
        canRunOnFailed: true, consoleParsers: [[parserName: 'Doxygen']],
        defaultEncoding: '', excludePattern: '', healthy: '',
        includePattern: '', messagesPattern: '', unHealthy: '',
        unstableNewAll: '0', useStableBuildAsReference: true])
}

if (helper.isOriginMaster(this)) {
    stage('Deploy (Docs)') {
        sh 'rsync -a --delete --stats build/docs/ ' +
            'web@doxygen.opengeosys.org:/www/doxygenogs'
    }
}
