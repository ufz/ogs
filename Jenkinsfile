#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.7') _

def builders = [:]
def helper = new ogs.helper()

timestamps {

builders['gcc'] = {
    node('docker') {
        dir('ogs') { checkoutWithTags() }
        load 'ogs/scripts/jenkins/gcc.groovy'
    }
}

builders['gcc-conan'] = {
    node('docker') {
        dir('ogs') { checkoutWithTags() }
        load 'ogs/scripts/jenkins/gcc-conan.groovy'
    }
}

builders['gcc-dynamic'] = {
    node('envinf1') {
        dir('ogs') { checkoutWithTags() }
        load 'ogs/scripts/jenkins/gcc-dynamic.groovy'
    }
}

builders['msvc'] = {
    node('win && conan') {
        dir('ogs') { checkoutWithTags() }
        load 'ogs/scripts/jenkins/msvc.groovy'
    }
}

builders['mac'] = {
    node('mac') {
        dir('ogs') { checkoutWithTags() }
        // load 'ogs/scripts/jenkins/mac.groovy'
    }
}

if (helper.isRelease(this)) {
    builders['msvc32'] = {
        node('win && conan') {
            dir('ogs') { checkoutWithTags() }
            load 'ogs/scripts/jenkins/msvc32.groovy'
        }
    }
}

builders['docs'] = {
    node('docker') {
        dir('ogs') { checkoutWithTags() }
        load 'ogs/scripts/jenkins/docs.groovy'
    }
}

parallel builders

def tag = ""
node('master') {
    checkoutWithTags()
    tag = helper.getTag()
    step([$class: 'LogParserPublisher',
        failBuildOnError: true,
        projectRulePath: "scripts/jenkins/all-log-parser.rules",
        showGraphs: true,
        unstableOnWarning: false,
        useProjectRule: true
    ])
}

if (currentBuild.result == "SUCCESS" || currentBuild.result == "UNSTABLE") {
    if (helper.isOriginMaster(this)) {
        build job: 'OGS-6/clang-sanitizer', wait: false
        if (tag != "") {
            keepBuild()
            currentBuild.displayName = tag
            helper.notification(msg: "Marked build for ${tag}.", script: this)

            node('mac') {
                dir('ogs') { checkoutWithTags() }
                load 'ogs/scripts/jenkins/docset.groovy'
            }
        }
    }
} else {
    if (helper.isOriginMaster(this)) {
        helper.notification(title: "${env.JOB_NAME} failed!", script: this,
            msg: "Build failed", url: "${env.BUILD_URL}/flowGraphTable/")
    }
}


} // end timestamps
