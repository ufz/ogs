def defaultCMakeOptions =
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System ' +
    '-DDOCS_GENERATE_LOGFILE=ON'

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

def image = docker.image('ogs6/gcc-latex:latest')
image.pull()
image.inside() {
    stage('Configure (Docs)') {
        configure.linux(cmakeOptions: defaultCMakeOptions, script: this)
    }

    stage('Generate (Docs)') {
        build.linux(script: this, target: 'doc')
    }
}

stage('Reports (Docs)') {
    publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: true,
        keepAll: true, reportDir: 'build/docs', reportFiles: 'index.html',
        reportName: 'Doxygen'])
    step([$class: 'WarningsPublisher', canResolveRelativePaths: false,
        messagesPattern:
            '.*DOT_GRAPH_MAX_NODES.,' +
            '.*potential recursive class relation.*',
        parserConfigurations: [[parserName: 'Doxygen', pattern:
        'build/DoxygenWarnings.log']], unstableTotalAll: '0'])
}

if (helper.isOriginMaster(this)) {
    stage('Deploy (Docs)') {
        sshagent(credentials: ['www-data_jenkins']) {
            sh 'rsync -a --delete --stats -e "ssh -o StrictHostKeyChecking=no"' +
                ' build/docs/ www-data@jenkins.opengeosys.org:'+
                '/var/www/doxygen.opengeosys.org'
        }
    }
}
