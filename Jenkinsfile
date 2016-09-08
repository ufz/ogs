#!/usr/bin/env groovy

node('master') {

    step([$class: 'GitHubSetCommitStatusBuilder', statusMessage: [content: 'Started Jenkins build']])
    checkout scm

    def builders = [:]
    builders['gcc'] = { load 'scripts/jenkins/gcc.groovy' }
    builders['msvc'] = { load 'scripts/jenkins/msvc.groovy' }
    builders['msvc32'] = { load 'scripts/jenkins/msvc32.groovy' }

    parallel builders

    step([$class: 'GitHubCommitStatusSetter'])
}

properties ([[$class: 'BuildDiscarderProperty', strategy: [$class: 'LogRotator', artifactDaysToKeepStr: '', artifactNumToKeepStr: '5', daysToKeepStr: '', numToKeepStr: '25']]])
