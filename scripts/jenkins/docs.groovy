node('docker') {
    def defaultCMakeOptions =
        '-DOGS_LIB_BOOST=System ' +
        '-DOGS_LIB_VTK=System'

    stage 'Checkout (Docs)'
    dir('ogs') { checkout scm }

    docker.image('ogs6/gcc-gui:latest').inside() {
        stage 'Configure (Docs)'
        configure.linux 'build', "${defaultCMakeOptions}"

        stage 'Generate (Docs)'
        build.linux 'build', 'doc'
    }

    stage 'Reports (Docs)'
    publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: false, keepAll: false,
        reportDir: 'build/docs', reportFiles: 'index.html', reportName: 'Doxygen'])
    step([$class: 'WarningsPublisher', canComputeNew: false,
        canResolveRelativePaths: false, consoleParsers: [[parserName: 'Doxygen']],
        defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '',
        messagesPattern: '', unHealthy: ''])

    if (helper.isOriginMaster()) {
        stage 'Deploy (Docs)'
        sh 'rsync -a --delete --stats build/docs/ ' +
            'web@doxygen.opengeosys.org:/www/doxygenogs'
    }
}
