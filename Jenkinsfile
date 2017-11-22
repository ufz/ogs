#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.9') _

pipeline {
  agent none
  options {
    ansiColor('xterm')
    timestamps()
  }
  stages {
    stage('Build') {
      parallel {
        // ************************** Docker ***********************************
        stage('Docker') {
          agent {
            docker {
              image 'ogs6/gcc-gui:latest'
              label 'docker'
              args '-v /home/jenkins/.ccache:/usr/src/.ccache'
              alwaysPull true
            }
          }
          environment {
            CONTENTFUL_ACCESS_TOKEN = credentials('CONTENTFUL_ACCESS_TOKEN')
            CONTENTFUL_OGS_SPACE_ID = credentials('CONTENTFUL_OGS_SPACE_ID')
          }
          steps {
            // Install web dependencies
            script {
              sh("""
                cd web
                yarn --ignore-engines --non-interactive
                node node_modules/node-sass/scripts/install.js
                npm rebuild node-sass
                sudo -H pip install -r requirements.txt
                """.stripIndent())

              configure { cmakeOptions =
                '-DOGS_CPU_ARCHITECTURE=generic ' +
                '-DOGS_WEB_BASE_URL=$JOB_URL"Web/" ' // TODO: or '-DOGS_WEB_BASE_URL=https://dev.opengeosys.org'
              }
              build { }
              build { target="tests" }
              build { target="ctest" }
              build { target="web" }
              build { target="doc" }
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
                dir('web/public') { stash(name: 'web') }
                dir('build/docs') { stash(name: 'doxygen') }
                script {
                  publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: true,
                    keepAll: true, reportDir: 'build/docs', reportFiles: 'index.html',
                    reportName: 'Doxygen'])
                  publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: true,
                    keepAll: true, reportDir: 'web/public', reportFiles: 'index.html',
                    reportName: 'Web'])
                  step([$class: 'WarningsPublisher', canResolveRelativePaths: false,
                    messagesPattern: """
                      .*DOT_GRAPH_MAX_NODES.
                      .*potential recursive class relation.*""",
                    parserConfigurations: [[parserName: 'Doxygen', pattern:
                    'build/DoxygenWarnings.log']], unstableTotalAll: '0'])
                }
                dir('build') { deleteDir() }
            }
          }
        }
        // ************************ Docker-Conan *******************************
        stage('Docker-Conan') {
          agent {
            docker {
              image 'ogs6/gcc-conan:latest'
              label 'docker'
              args '-v /home/jenkins/.ccache:/usr/src/.ccache'
              alwaysPull true
            }
          }
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DOGS_USE_CONAN=ON ' +
                  '-DOGS_CONAN_BUILD=never ' +
                  '-DOGS_CPU_ARCHITECTURE=generic ' +
                  '-DOGS_PACKAGE_DEPENDENCIES=ON '
              }
              build { }
              build { target="tests" }
              build { target="ctest" }
              configure {
                cmakeOptions =
                  '-DOGS_BUILD_CLI=OFF ' +
                  '-DOGS_USE_PCH=OFF ' +     // see #1992
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_TESTS=OFF ' +
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
                archiveArtifacts 'build/*.tar.gz,build/conaninfo.txt'
                dir('build') { deleteDir() }
            }
          }
        }
        // ************************** envinf1 **********************************
        stage('Envinf1 (serial)') {
          agent { label "envinf1"}
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_METIS=ON ' +
                  '-DBUILD_SHARED_LIBS=ON '
                env = 'envinf1/cli.sh'
              }
              build {
                env = 'envinf1/cli.sh'
                cmd_args = '-l 30'
              }
              build {
                env = 'envinf1/cli.sh'
                target = 'tests'
              }
              build {
                env = 'envinf1/cli.sh'
                target = 'ctest'
              }
            }
          }
          post {
            always {
              publishReports { }
              dir('build') { deleteDir() }
            }
          }
        }
        stage('Envinf1 (parallel)') {
          agent { label "envinf1"}
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_METIS=ON ' +
                  '-DBUILD_SHARED_LIBS=ON ' +
                  '-DOGS_USE_PETSC=ON '
                env = 'envinf1/petsc.sh'
                generator = 'Unix Makefiles'
              }
              build {
                env = 'envinf1/petsc.sh'
                cmd_args = '-l 30'
              }
              build {
                env = 'envinf1/petsc.sh'
                target = 'tests'
              }
              build {
                env = 'envinf1/petsc.sh'
                target = 'ctest'
              }
            }
          }
          post {
            always {
              publishReports { }
              dir('build') { deleteDir() }
            }
          }
        }
        // ************************** Windows **********************************
        stage('Win') {
          agent {label 'win && conan' }
          environment {
            MSVC_NUMBER = '15'
            MSVC_VERSION = '2017'
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
      } // end parallel
    } // end stage Build
    // *************************** Log Parser **********************************
    stage('Log Parser') {
      agent any
      steps {
        script {
          checkout scm
          step([$class: 'LogParserPublisher',
              failBuildOnError: true,
              projectRulePath: "scripts/jenkins/all-log-parser.rules",
              showGraphs: true,
              unstableOnWarning: false,
              useProjectRule: true
          ])
        }
      }
    }
    stage('Master') {
      when { environment name: 'JOB_NAME', value: 'ufz/ogs/master' }
      parallel {
        // ************************* Deploy Web ********************************
        stage('Deploy Web') {
          agent any
          steps {
            dir('web') { unstash 'web' }
            dir('doxygen') { unstash 'doxygen' }
            script {
              sshagent(credentials: ['www-data_jenkins']) {
                sh 'rsync -a --delete --stats -e "ssh -o UserKnownHostsFile=' +
                   'ogs/scripts/jenkins/known_hosts" . ' +
                   'www-data@jenkins.opengeosys.org:/var/www/dev.opengeosys.org'
                sh 'rsync -a --delete --stats -e "ssh -o UserKnownHostsFile=' +
                   'ogs/scripts/jenkins/known_hosts" . ' +
                   'www-data@jenkins.opengeosys.org:/var/www/doxygen.opengeosys.org'
              }
            }
          }
        }
        // *********************** Deploy envinf1 ******************************
        stage('Deploy envinf1') {
          agent { label "envinf1"}
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_METIS=ON ' +
                  '-DBUILD_SHARED_LIBS=ON ' +
                  '-DCMAKE_INSTALL_PREFIX=/global/apps/ogs/head/standard ' +
                  '-DOGS_MODULEFILE=/global/apps/modulefiles/ogs/head/standard ' +
                  '-DOGS_CPU_ARCHITECTURE=core-avx-i '
                env = 'envinf1/cli.sh'
              }
              build {
                env = 'envinf1/cli.sh'
                target = 'install'
                cmd_args = '-l 30'
              }
            }
          }
          post {
            always {
              dir('build') { deleteDir() }
            }
          }
        }
        // ******************** Deploy envinf1 PETSc ***************************
        stage('Deploy envinf1 PETSc') {
          agent { label "envinf1"}
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DOGS_USE_PETSC=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_METIS=ON ' +
                  '-DBUILD_SHARED_LIBS=ON ' +
                  '-DCMAKE_INSTALL_PREFIX=/global/apps/ogs/head/petsc ' +
                  '-DOGS_MODULEFILE=/global/apps/modulefiles/ogs/head/petsc ' +
                  '-DOGS_CPU_ARCHITECTURE=core-avx-i '
                env = 'envinf1/petsc.sh'
              }
              build {
                env = 'envinf1/petsc.sh'
                target = 'install'
                cmd_args = '-l 30'
              }
            }
          }
          post {
            always {
              dir('build') { deleteDir() }
            }
          }
        }
        // ************************** Sanitizer ********************************
        stage('Sanitizer') {
          agent {
            docker {
              image 'ogs6/clang-base:latest'
              label 'docker'
              args '-v /home/jenkins/.ccache:/usr/src/.ccache'
              alwaysPull true
            }
          }
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DOGS_ADDRESS_SANITIZER=ON ' +
                  '-DOGS_UNDEFINED_BEHAVIOR_SANITIZER=ON ' +
                  '-DOGS_BUILD_UTILS=ON '
              }
              try {
                build {
                  target = 'test'
                  cmd = 'UBSAN_OPTIONS=print_stacktrace=1 make -j $(nproc)'
                }
              }
              catch(err) { echo "Clang sanitizer for unit tests failed!" }

              try {
                build {
                  target = 'ctest'
                  cmd = 'UBSAN_OPTIONS=print_stacktrace=1 make -j $(nproc)'
                }
              }
              catch(err) { echo "Clang sanitizer for end-to-end tests failed!" }
            }
          }
          post {
            always {
              dir('build') { deleteDir() }
            }
          }
        }
        // *************************** Post ************************************
        stage('Post') {
          agent any
          steps {
            script {
              def helper = new ogs.helper()
              checkout scm
              def tag = helper.getTag()
              if (tag != "") {
                keepBuild()
                currentBuild.displayName = tag
                helper.notification(msg: "Marked build for ${tag}.", script: this)
              }
            }
          }
        }
      } // end parallel
    } // end stage master
  }
}
