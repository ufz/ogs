#!/usr/bin/env groovy

node('master') {
    step([$class: 'GitHubSetCommitStatusBuilder', statusMessage: [content: 'Started Jenkins build']])
    checkout scm

    def builders = [:]
    builders['gcc'] = { load 'scripts/jenkins/gcc.groovy' }
    builders['msvc'] = { load 'scripts/jenkins/msvc.groovy' }
    builders['mac'] = { load 'scripts/jenkins/mac.groovy' }

    helper = load 'scripts/jenkins/lib/helper.groovy'
    if (helper.isRelease()) {
        builders['msvc32'] = { load 'scripts/jenkins/msvc32.groovy' }
    }

    parallel builders

    step([$class: 'GitHubCommitStatusSetter'])

    if (helper.isRelease()) { build job: 'OGS-6/Deploy', wait: false }
}

properties ([[$class: 'BuildDiscarderProperty', strategy: [$class: 'LogRotator', artifactDaysToKeepStr: '', artifactNumToKeepStr: '5', daysToKeepStr: '', numToKeepStr: '25']]])
