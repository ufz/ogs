#!/usr/bin/env groovy
def builders = [:]
def helper = new ogs.helper()

timestamps {

builders['gcc'] = {
    node('docker') {
        dir('ogs') { checkout scm }
        load 'ogs/scripts/jenkins/gcc.groovy'
    }
}

builders['gcc-dynamic'] = {
    node('docker') {
        dir('ogs') { checkout scm }
        load 'ogs/scripts/jenkins/gcc-dynamic.groovy'
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
        def tag = ""
        node('master') {
            checkout scm
            tag = helper.getTag()
        }
        if (tag != "") {
            keepBuild()
            currentBuild.displayName = tag

            node('mac') {
                dir('ogs') { checkout scm }
                load 'ogs/scripts/jenkins/docset.groovy'
            }
        }
    }
}

} // end timestamps
