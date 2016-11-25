def defaultCMakeOptions =
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System ' +
    '-DDOCS_GENERATE_LOGFILE=ON'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

def image = docker.image('ogs6/gcc-base:latest')
image.pull()
image.inside() {
    stage('Configure (Docs)') {
        configure.linux 'build', "${defaultCMakeOptions}"
    }

    stage('Generate (Docs)') {
        build.linux this, 'build', 'doc'
    }
}

stage('Reports (Docs)') {
    publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: false,
        keepAll: false, reportDir: 'build/docs', reportFiles: 'index.html',
        reportName: 'Doxygen'])
    step([$class: 'WarningsPublisher', canResolveRelativePaths: false,
        parserConfigurations: [[parserName: 'Doxygen', pattern:
        'build/DoxygenWarnings.log']], unstableNewAll: '0'])
}

if (helper.isOriginMaster(this)) {
    stage('Deploy (Docs)') {
        sh 'rsync -a --delete --stats build/docs/ ' +
            'web@doxygen.opengeosys.org:/www/doxygenogs'
    }
}
