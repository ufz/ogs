defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
defaultCMakeOptions = '-DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System'

node('docker') {
    stage 'Checkout'
    checkout scm

    stage 'Build'
    docker.image('ogs6/gcc-base:latest').inside(defaultDockerArgs) {
        build 'build', '-DOGS_COVERAGE=ON',
            'testrunner_coverage_cobertura ctest_coverage_cobertura'
    }

    archive 'build/*.xml'
    // Report is published in a free-style child job
}

def build(buildDir, cmakeOptions, target) {
    sh "rm -rf ${buildDir} && mkdir ${buildDir}"

    stage 'Configure'
    sh "cd ${buildDir} && cmake ../ogs ${defaultCMakeOptions} ${cmakeOptions}"

    stage 'Build'
    sh "cd ${buildDir} && make -j \$(nproc) ${target}"
}

properties([[$class: 'BuildDiscarderProperty', strategy: [$class: 'LogRotator',
    artifactDaysToKeepStr: '', artifactNumToKeepStr: '1', daysToKeepStr: '', numToKeepStr: '5']]])
