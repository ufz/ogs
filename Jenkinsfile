#!/usr/bin/env groovy
@Library('jenkins-pipeline@master') _

pipeline {
  agent none
  options {
    ansiColor('xterm')
    timestamps()
    // alwaysDoDockerPull()
  }
  stages {
    stage('Build') {
      parallel {
        stage('Docker') {
          agent {
            docker {
              image 'ogs6/gcc-gui:latest'
              label 'docker'
              args '-v /home/jenkins/.ccache:/usr/src/.ccache'
            }
          }
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DCMAKE_BUILD_TYPE=Release ' +
                  '-DOGS_CPU_ARCHITECTURE=generic ' +
                  '-DOGS_WEB_BASE_URL=$JOB_URL"Web/" '
              }
              build { }
              build { target="tests" }
              build { target="ctest" }
            }
          }
        }
        stage('Win') {
          agent {
            label "win && conan"
          }
          environment {
            MSVC_NUMBER = '15'
            MSVC_VERSION = "2017"
          }
          steps {
            script {
              // CLI
              configure {
                cmakeOptions =
                  '-DOGS_USE_CONAN=ON ' +
                  '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON '
              }
              build { }
              build { target="tests" }
              build { target="ctest" }
              // GUI
              configure {
                cmakeOptions =
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_TESTS=OFF ' +
                  '-DOGS_BUILD_SWMM=ON ' +
                  '-DOGS_BUILD_METIS=ON '
                  keepDir = true
              }
              build { }
            }
          }
          post {
            always {
              publishReports { }
            }
            failure {
                dir('build') { deleteDir() }
            }
            success {
                archiveArtifacts 'build/*.zip,build/conaninfo.txt'
                dir('build') { deleteDir() }
            }
          }
        }
      }
    }
  }
}
