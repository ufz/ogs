#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.4') _

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

if (currentBuild.result == "SUCCESS" || currentBuild.result == "UNSTABLE") {
    if (helper.isOriginMaster(this)) {
        build job: 'OGS-6/clang-sanitizer', wait: false
        def tag = ""
        node('master') {
            checkoutWithTags()
            tag = helper.getTag()
        }
        if (tag != "") {
            keepBuild()
            currentBuild.displayName = tag

            node('mac') {
                dir('ogs') { checkoutWithTags() }
                load 'ogs/scripts/jenkins/docset.groovy'
            }
        }
    }
}

} // end timestamps
