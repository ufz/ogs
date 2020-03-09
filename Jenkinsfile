#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.22') _

def stage_required = [build: false, full: false]
def build_shared = 'ON'

pipeline {
  agent none
  options {
    ansiColor('xterm')
    timestamps()
    buildDiscarder(logRotator(numToKeepStr: '30', artifactNumToKeepStr: '10'))
    timeout(time: 6, unit: 'HOURS')
  }
  parameters {
    booleanParam(name: 'docker_conan', defaultValue: true)
    booleanParam(name: 'docker_conan_debug', defaultValue: true)
    booleanParam(name: 'docker_conan_gui', defaultValue: true)
    booleanParam(name: 'eve_serial', defaultValue: true)
    booleanParam(name: 'eve_parallel', defaultValue: true)
    booleanParam(name: 'win', defaultValue: true)
    booleanParam(name: 'mac', defaultValue: true)
    booleanParam(name: 'clang_analyzer', defaultValue: true)
    booleanParam(name: 'master_jobs', defaultValue: true)
  }
  stages {
    stage('Build') {
      parallel {
        // ************************** Windows **********************************
        stage('Win') {
          agent {label 'win && conan' }
          environment {
            MSVC_NUMBER = '15'
            MSVC_VERSION = '2017'
            OMP_NUM_THREADS = '1'
            CC = 'clcache'
            CXX = 'clcache'
          }
          steps {
            script {
              def num_threads = env.NUM_THREADS
              bat 'git submodule sync && git submodule update'
              bat 'conan remove --locks'
              configure { // CLI + GUI
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=OFF " +
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_CONAN_BUILD=missing ' +
                  '-DOGS_BUILD_SWMM=ON ' +
                  '-DOGS_USE_NETCDF=ON '
              }
              build {
                target="MeshLib"
                log="build1.log"
                cmd_args="-l ${num_threads}"
              }
            }
          }
          post {
            always {
              xunit([
                CTest(pattern: 'build/Testing/**/*.xml'),
                GoogleTest(pattern: 'build/Tests/testrunner.xml')
              ])
              recordIssues enabledForFailure: true, filters: [
                excludeFile('.*\\.conan.*'), excludeFile('.*ThirdParty.*'),
                excludeFile('.*thread.hpp')],
                tools: [msBuild(name: 'MSVC', pattern: 'build/build*.log')],
                qualityGates: [[threshold: 10, type: 'TOTAL', unstable: true]]
            }
            success {
              archiveArtifacts 'build/*.zip,build/conaninfo.txt'
            }
          }
        }
      } // end parallel
    } // end stage Build
  }
}
