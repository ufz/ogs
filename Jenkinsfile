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
     // *************************** Git Check **********************************
    stage('Pre-checks') {
      agent { label "master"}
      steps {
        sh "pre-commit install"
        sh "pre-commit run --all-files"
        sh "git config core.whitespace -blank-at-eof"
        sh "git diff --check `git merge-base origin/master HEAD` HEAD -- . ':!*.md' ':!*.pandoc' ':!*.asc' ':!*.dat' ':!*.ts'"
        dir('scripts/jenkins') { stash(name: 'known_hosts', includes: 'known_hosts') }
        ciSkip action: 'check' // Check for [ci skip] or [web] commit message.

        // ********* Check changesets for conditional stage execution **********
        script {
          def causes = currentBuild.getBuildCauses()
          for(cause in causes) {
            if (cause.class.toString().contains("UserIdCause")) {
              echo "Doing full build because job was started by user."
              stage_required.full = true
              env.CI_SKIP = "false"
              return true
            }
          }

          if (env.JOB_NAME == 'ufz/ogs/master') {
            build_shared = 'OFF'
          }
          if (currentBuild.number == 1 || buildingTag()) {
            stage_required.full = true
            return true
          }
          def changeLogSets = currentBuild.changeSets
          for (int i = 0; i < changeLogSets.size(); i++) {
            def entries = changeLogSets[i].items
            for (int j = 0; j < entries.length; j++) {
              def paths = new ArrayList(entries[j].affectedPaths)
              for (int p = 0; p < paths.size(); p++) {
                def path = paths[p]
                if (path.matches("Jenkinsfile")) {
                  stage_required.full = true
                  echo "Doing full build."
                  return true
                }
                if (path.matches("^(CMakeLists.txt|scripts|Applications|BaseLib|FileIO|GeoLib|MaterialLib|MathLib|MeshGeoToolsLib|MeshLib|NumLib|ProcessLib|SimpleTests|Tests).*")
                  && !stage_required.build) {
                  stage_required.build = true
                  echo "Doing regular build."
                }
              }
            }
          }
          if(!(stage_required.build || stage_required.full)) {
            currentBuild.result='NOT_BUILT'
          }
        }
      }
      // Mark build as NOT_BUILT when [ci skip] commit message was found
      post { always { ciSkip action: 'postProcess' } }
    }
    stage('Build') {
      parallel {
        // **************************** eve ************************************
        stage('Frontend2 (serial)') {
          when {
            beforeAgent true
            expression { return params.eve_serial && (stage_required.build || stage_required.full) }
          }
          agent { label "frontend2"}
          environment {
            OMP_NUM_THREADS = '1'
          }
          steps {
            script {
              sh 'git submodule sync'
              configure {
                cmakeOptions =
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DBUILD_SHARED_LIBS=ON ' +
                  '-DOGS_USE_PYTHON=ON ' +
                  '-DCMAKE_INSTALL_PREFIX=/global/apps/ogs/head/standard ' +
                  '-DOGS_MODULEFILE=/global/apps/modulefiles/ogs/head/standard '
                env = 'eve/cli.sh'
              }
              build {
                env = 'eve/cli.sh'
                cmd_args = '-l 30'
              }
              build {
                env = 'eve/cli.sh'
                target = 'tests'
              }
              build {
                env = 'eve/cli.sh'
                target = 'ctest'
              }
            }
          }
          post {
            always {
              xunit([
                CTest(pattern: 'build/Testing/**/*.xml'),
                GoogleTest(pattern: 'build/Tests/testrunner.xml')
              ])
            }
            success {
              script {
                if (env.JOB_NAME == 'ufz/ogs/master') {
                  sh 'rm -rf /global/apps/ogs/head/standard'
                  build {
                    env = 'eve/cli.sh'
                    target = 'install'
                  }
                }
              }
            }
          }
        }
      } // end parallel
    } // end stage Build
  }
}
