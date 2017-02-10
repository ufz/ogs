def defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
def defaultCMakeOptions =
    '-DCMAKE_BUILD_TYPE=Release ' +
    '-DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_EIGEN=Local ' +
    '-DOGS_LIB_VTK=System '

def guiCMakeOptions =
    '-DOGS_BUILD_CLI=OFF ' +
    '-DOGS_BUILD_GUI=ON ' +
    '-DOGS_BUILD_UTILS=ON ' +
    '-DOGS_BUILD_TESTS=OFF ' +
    '-DOGS_BUILD_METIS=ON '

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

def webCMakeOptions = '-DOGS_WEB_BASE_URL=$JOB_URL"Web/" '
if (helper.isOriginMaster(this))
    webCMakeOptions = '-DOGS_WEB_BASE_URL=https://dev.opengeosys.org'

def image = docker.image('ogs6/gcc-gui:latest')
image.pull()
image.inside(defaultDockerArgs) {
    sh 'cd ogs && git lfs pull'
    stage('Install prerequisites Web') {
        sh 'cd ogs/web && rm -rf node_modules && yarn && sudo -H pip install -r requirements.txt'
    }
    stage('Configure (Linux-Docker)') {
        configure.linux(cmakeOptions: defaultCMakeOptions + webCMakeOptions, script: this)
    }

    stage('CLI (Linux-Docker)') {
        build.linux(script: this)
    }

    stage('Test (Linux-Docker)') {
        build.linux(script: this, target: 'tests ctest')
    }

    stage('Web (Linux-Docker)') {
        withCredentials([string(
            credentialsId: 'CONTENTFUL_ACCESS_TOKEN',
            variable: 'CONTENTFUL_ACCESS_TOKEN'), string(
            credentialsId: 'CONTENTFUL_OGS_SPACE_ID',
            variable: 'CONTENTFUL_OGS_SPACE_ID')]) {

            build.linux(script: this, target: 'web')

            if (helper.isOriginMaster(this)) {
                sshagent(credentials: ['www-data_jenkins']) {
                    sh 'rsync -a --delete --stats -e "ssh -o StrictHostKeyChecking=no"' +
                        ' ogs/web/public/ www-data@jenkins.opengeosys.org:'+
                        '/var/www/dev.opengeosys.org'
                }
            } else {
                publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: false,
                    keepAll: false, reportDir: 'ogs/web/public', reportFiles: 'index.html',
                    reportName: 'Web'])
            }
        }
    }

    stage('Data Explorer (Linux-Docker)') {
        configure.linux(
            cmakeOptions: defaultCMakeOptions + guiCMakeOptions,
            keepDir: true,
            script: this
        )
        build.linux(script: this)
    }
}

stage('Post (Linux-Docker)') {
    post.publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/clang-log-parser.rules'
    post.cleanup()
}
