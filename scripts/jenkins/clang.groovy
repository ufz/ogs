defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'

node('docker') {
    stage 'Checkout'
    dir('ogs') { checkout scm }

    // Multiple configurations are build in parallel
    parallel linux: {
        docker.image('ogs6/clang-base:latest').inside(defaultDockerArgs) {
            catchError {
                build 'build', '-DOGS_ADDRESS_SANITIZER=ON ' +
                    '-DOGS_UNDEFINED_BEHAVIOR_SANITIZER=ON', ''

                stage 'Unit tests'
                sh '''cd build
                      UBSAN_OPTIONS=print_stacktrace=1 make tests'''

                stage 'End-to-end tests'
                sh '''cd build
                      UBSAN_OPTIONS=print_stacktrace=1 make ctest'''
            }
        }
        step([$class: 'LogParserPublisher', failBuildOnError: true, unstableOnWarning: true,
            projectRulePath: 'ogs/scripts/jenkins/clang-log-parser.rules', useProjectRule: true])
    }

    step([$class: 'JUnitResultArchiver',
        testResults: 'build/Tests/testrunner.xml'])

    archive 'build*/*.tar.gz'
} // end node


def build(buildDir, cmakeOptions, target) {
    sh "rm -rf ${buildDir} && mkdir ${buildDir}"

    stage 'Configure'
    sh "cd ${buildDir} && cmake ../ogs ${cmakeOptions}"

    stage 'Build'
    sh "cd ${buildDir} && make -j \$(nproc) ${target}"
}

properties([[$class: 'BuildDiscarderProperty', strategy: [$class: 'LogRotator',
    artifactDaysToKeepStr: '', artifactNumToKeepStr: '5', daysToKeepStr: '', numToKeepStr: '25']]])
