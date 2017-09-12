def defaultCMakeOptions =
    '-DCMAKE_BUILD_TYPE=Release ' +
    '-DOGS_BUILD_UTILS=ON ' +
    '-DOGS_BUILD_METIS=ON '

def configure = new ogs.configure()
def build = new ogs.build()
def post = new ogs.post()
def helper = new ogs.helper()

stage('Configure (envinf1)') {
    installPrefix = "/global/apps/ogs/"
    installPrefix += "head"
    modulePrefix = "/global/apps/modulefiles/ogs/"
    modulePrefix += "head"

    configure.linux(cmakeOptions: defaultCMakeOptions + '-DBUILD_SHARED_LIBS=ON ',
        env: 'envinf1/cli.sh', script: this)
    configure.linux(cmakeOptions: defaultCMakeOptions + '-DBUILD_SHARED_LIBS=ON ' +
        '-DOGS_USE_PETSC=ON ',
        dir: 'build-petsc', env: 'envinf1/petsc.sh', script: this)
}

stage('CLI (envinf1)') {
    parallel serial: {
        build.linux(env: 'envinf1/cli.sh', script: this)
    }, petsc: {
        build.linux(dir: 'build-petsc', env: 'envinf1/petsc.sh', script: this)
    }
}

stage('Test (envinf1)') {
    parallel serial: {
        build.linux(env: 'envinf1/cli.sh', script: this, target: 'tests ctest')
    }, petsc: {
        build.linux(dir: 'build-petsc', env: 'envinf1/petsc.sh', script: this,
            target: 'tests ctest')
    }
}

if (helper.isOriginMaster(this)) {
    stage('Deploy (envinf1)') {
        parallel serial: {
            configure.linux(cmakeOptions: defaultCMakeOptions +
                "-DCMAKE_INSTALL_PREFIX=${installPrefix}/standard " +
                "-DOGS_MODULEFILE=${modulePrefix}/standard " +
                "-DOGS_CPU_ARCHITECTURE=core-avx-i ",
                dir: 'build-static', env: 'envinf1/cli.sh', script: this)
            build.linux(dir: 'build-static', env: 'envinf1/cli.sh',
                script: this, target: 'install')

        }, petsc: {
            configure.linux(cmakeOptions: defaultCMakeOptions + '-DOGS_USE_PETSC=ON ' +
                "-DCMAKE_INSTALL_PREFIX=${installPrefix}/petsc " +
                "-DOGS_MODULEFILE=${modulePrefix}/petsc " +
                "-DOGS_CPU_ARCHITECTURE=core-avx-i ",
                dir: 'build-static-petsc', env: 'envinf1/petsc.sh', script: this)
            build.linux(dir: 'build-static-petsc', env: 'envinf1/petsc.sh',
                script: this, target: 'install')
        }
    }
}

stage('Post (envinf1)') {
    post.publishTestReports 'build*/Testing/**/*.xml', 'build*/Tests/testrunner.xml'
    post.cleanup()
}
