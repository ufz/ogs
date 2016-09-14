#!/usr/bin/env groovy

node('master') {

    step([$class: 'GitHubSetCommitStatusBuilder', statusMessage: [content: 'Started Jenkins build']])
    checkout scm

    def builders = [:]
    builders['gcc'] = { load 'scripts/jenkins/gcc.groovy' }
    builders['msvc'] = { load 'scripts/jenkins/msvc.groovy' }
    builders['msvc32'] = { load 'scripts/jenkins/msvc32.groovy' }
    builders['mac'] = { load 'scripts/jenkins/mac.groovy' }

    parallel builders

    step([$class: 'GitHubCommitStatusSetter'])

    if (env.BRANCH_NAME == 'master') {
        build job: 'OGS-6/Deploy', wait: false
    }
}

properties ([[$class: 'BuildDiscarderProperty', strategy: [$class: 'LogRotator', artifactDaysToKeepStr: '', artifactNumToKeepStr: '5', daysToKeepStr: '', numToKeepStr: '25']]])
