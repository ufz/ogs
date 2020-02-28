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
        // ************************ Docker-Conan *******************************
        stage('Docker-Conan') {
          when {
            beforeAgent true
            expression { return params.docker_conan && (stage_required.build || stage_required.full) }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.gcc.full'
              dir 'scripts/docker'
              label 'docker'
              args '-v /home/jenkins/cache/ccache:/opt/ccache -v /home/jenkins/cache/conan/.conan:/opt/conan/.conan'
              additionalBuildArgs '--pull'
            }
          }
          environment {
            OMP_NUM_THREADS = '1'
            LD_LIBRARY_PATH = "$WORKSPACE/build/lib"
          }
          steps {
            script {
              sh 'git submodule sync && git submodule update'
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=OFF " +
                  '-DOGS_CPU_ARCHITECTURE=generic ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_CONAN_BUILD=missing ' +
                  '-DOGS_USE_CVODE=ON ' +
                  '-DOGS_USE_MFRONT=ON ' +
                  '-DOGS_USE_PYTHON=ON '
              }
              build {
                target="package"
                log="build1.log"
              }
              build { target="tests" }
              build { target="ctest" }
              build { target="doc" }
            }
            stash name: 'bats-files', includes: 'build/*.bats'
          }
          post {
            always {
              xunit([
                // Testing/-folder is a CTest convention
                CTest(pattern: 'build/Testing/**/*.xml'),
                GoogleTest(pattern: 'build/Tests/testrunner.xml')
              ])
              recordIssues enabledForFailure: true, filters: [
                excludeFile('.*qrc_icons\\.cpp.*'), excludeFile('.*QVTKWidget.*'),
                excludeFile('.*\\.conan.*/TFEL/.*'),
                excludeFile('.*ThirdParty/MGIS/src.*'),
                excludeFile('.*MaterialLib/SolidModels/MFront/.*\\.mfront'),
                excludeFile('.*MaterialLib/SolidModels/MFront/.*\\.hxx'),
                excludeMessage('.*tmpnam.*')],
                tools: [gcc4(name: 'GCC', pattern: 'build/build*.log')],
                unstableTotalAll: 1
              recordIssues enabledForFailure: true, filters: [
                  excludeFile('-'), excludeFile('.*Functional\\.h'),
                  excludeFile('.*gmock-.*\\.h'), excludeFile('.*gtest-.*\\.h')
                ],
                // Doxygen is handled by gcc4 parser as well
                tools: [gcc4(name: 'Doxygen', id: 'doxygen',
                             pattern: 'build/DoxygenWarnings.log')],
                unstableTotalAll: 1
            }
            success {
              publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: true,
                  keepAll: true, reportDir: 'build/docs', reportFiles: 'index.html',
                  reportName: 'Doxygen'])
              archiveArtifacts 'build/*.tar.gz,build/conaninfo.txt,build/*.bats'
              script {
                if (env.JOB_NAME == 'ufz/ogs/master') {
                  dir 'build' { stash name: 'bats-files', includes: '*.bats' }
                  // Deploy Doxygen
                  unstash 'known_hosts'
                  sshagent(credentials: ['www-data_jenkins']) {
                    sh 'rsync -a --delete --stats -e "ssh -o UserKnownHostsFile=' +
                       'known_hosts" build/docs/. ' +
                       'www-data@jenkins.opengeosys.org:/var/www/doxygen.opengeosys.org'
                  }
                }
              }
            }
          }
        }
        // ******************** Container Maker Images *************************
        stage('Container Maker Images') {
          when {
            beforeAgent true
            allOf {
              expression { return params.master_jobs }
              environment name: 'JOB_NAME', value: 'ufz/ogs/master'
            }
          }
          agent { label 'docker'}
          steps {
            script {
              sh '''git submodule update --init ThirdParty/container-maker
                virtualenv .venv
                source .venv/bin/activate
                pip install -r ThirdParty/container-maker/requirements.txt
                export PYTHONPATH="${PYTHONPATH}:${PWD}/ThirdParty/container-maker"
                python ThirdParty/container-maker/ogscm/cli.py -B -C -R \
                  -j $NUM_THREADS --ogs . --pm system --cvode \
                  --ompi off
                cp _out_images/ogs*serial*.sif ogs-serial.sif
              '''.stripIndent()
            }
            stash name: 'serial-sif', includes: 'ogs-serial.sif'
          }
          post {
            success {
              archiveArtifacts('_out/images/*.sif')
              dir('_out') { deleteDir() } // Cleanup
            }
          }
        }
      } // end parallel
    } // end stage Build
    stage('Post') {
      parallel {
        // ********************* External benchmarks ***************************
        stage('External Benchmarks') {
          agent { label "envinf1" }
          when {
            beforeAgent true
            environment name: 'JOB_NAME', value: 'ufz/ogs/master'
          }
          environment {
            SIF = "ogs-serial.sif"
            SRC = "$WORKSPACE/ogs-data"
            OUT = "$WORKSPACE/out"
          }
          stages {
            stage('Update ogs-data') {
              steps {
                script {
                  dir('ogs') { checkout scm }
                  dir('ogs-data') {
                    checkout(changelog: false, poll: false, scm: [$class: 'GitSCM',
                      extensions: [[$class: 'CloneOption', shallow: true]],
                      userRemoteConfigs: [[
                        credentialsId: '2719b702-1298-4e87-8464-5dfc62fbd923',
                        url: 'https://github.com/ufz/ogs-data']]])
                    sh 'rsync -av --delete --exclude .git/ ../ogs/Tests/Data/ .'
                    unstash 'bats-files'
                    sh 'mv build/*.bats . && rm -r build'
                    sh "git add --all . && git diff --quiet && git diff --staged --quiet || git commit -am 'Update'"
                    withCredentials([usernamePassword(
                      credentialsId: '2719b702-1298-4e87-8464-5dfc62fbd923',
                      passwordVariable: 'GIT_PASSWORD', usernameVariable: 'GIT_USERNAME')]) {
                      sh 'git push https://${GIT_USERNAME}:${GIT_PASSWORD}@github.com/ufz/ogs-data HEAD:master'
                    }
                  }
                }
              }
            }
            stage('Container benchmarks') {
              steps {
                unstash 'serial-sif'
                sh 'bats ogs-data/benchmarks.bats'
              }
            }
          }
        }
      } // end parallel
    } // end stage Post
  }
}
