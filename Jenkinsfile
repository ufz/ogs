#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.8') _

def builders = [:]
def helper = new ogs.helper()
def tag = ""

timestamps {

builders['gcc'] = {
    node('docker') {
        dir('ogs') {
            checkout scm
            tag = helper.getTag()
        }
        load 'ogs/scripts/jenkins/gcc.groovy'
    }
}

builders['gcc-conan'] = {
    node('docker') {
        dir('ogs') { checkout scm }
        load 'ogs/scripts/jenkins/gcc-conan.groovy'
    }
}

builders['gcc-dynamic'] = {
    node('envinf1') {
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

// builders['mac'] = {
    // node('mac') {
        // dir('ogs') { checkout scm }
        // load 'ogs/scripts/jenkins/mac.groovy'
    // }
// }

if (helper.isRelease(this)) {
    builders['msvc32'] = {
        node('win && conan') {
            dir('ogs') { checkout scm }
            load 'ogs/scripts/jenkins/msvc32.groovy'
        }
    }
}

if (helper.isOriginMaster(this)) {
    builders['coverage'] = {
        node('docker') {
            dir('ogs') { checkout scm }
            load 'ogs/scripts/jenkins/coverage.groovy'
        }
    }
}

try {
    parallel builders
}
catch (err) {
    currentBuild.result = 'FAILURE'
    if (helper.isOriginMaster(this)) {
        helper.notification(title: "${env.JOB_NAME} failed!", script: this,
            msg: "Build failed", url: "${env.BUILD_URL}/flowGraphTable/")
    }
}

node('master') {
    checkout scm
    step([$class: 'LogParserPublisher',
        failBuildOnError: true,
        projectRulePath: "scripts/jenkins/all-log-parser.rules",
        showGraphs: true,
        unstableOnWarning: false,
        useProjectRule: true
    ])
}

if (helper.isOriginMaster(this)) {
    if (currentBuild.result == "SUCCESS" || currentBuild.result == "UNSTABLE") {
        build job: 'OGS-6/clang-sanitizer', wait: false
        if (tag != "") {
            keepBuild()
            currentBuild.displayName = tag
            helper.notification(msg: "Marked build for ${tag}.", script: this)

            node('mac') {
                dir('ogs') { checkout scm }
                load 'ogs/scripts/jenkins/docset.groovy'
            }
        }
    }
}


} // end timestamps
