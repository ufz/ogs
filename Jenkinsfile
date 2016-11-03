#!/usr/bin/env groovy
def builders = [:]
def helper = new ogs.helper()

builders['gcc'] = {
    node('docker') {
        dir('ogs') { checkout scm }
        load 'ogs/scripts/jenkins/gcc.groovy'
    }
}

builders['msvc'] = {
    node('win && conan') {
        dir('ogs') { checkout scm }
        load 'ogs/scripts/jenkins/msvc.groovy'
    }
}

builders['mac'] = {
    node('mac') {
        dir('ogs') { checkout scm }
        load 'ogs/scripts/jenkins/mac.groovy'
    }
}

if (helper.isRelease(this)) {
    builders['msvc32'] = {
        node('win && conan') {
            dir('ogs') { checkout scm }
            load 'ogs/scripts/jenkins/msvc32.groovy'
        }
    }
}
if (helper.isOriginMaster(this)) {
    builders['docs'] = {
        node('docker') {
            dir('ogs') { checkout scm }
            load 'ogs/scripts/jenkins/docs.groovy'
        }
    }
}

parallel builders

node { step([$class: 'GitHubCommitStatusSetter']) }

if (currentBuild.result == "SUCCESS" || currentBuild.result == "UNSTABLE") {
    if (helper.isOriginMaster(this)) {
        build job: 'OGS-6/clang-sanitizer', wait: false
        node('master') {
            checkout scm
            def tag = helper.getTag()
            if (tag != "") {
                keepBuild()
                currentBuild.displayName = tag
            }
        }
    }
}
