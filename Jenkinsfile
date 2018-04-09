#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.9') _

def stage_required = [web: false, build: false, data: false, full: false]

pipeline {
  agent none
  options {
    ansiColor('xterm')
    timestamps()
  }
  stages {
     // *************************** Git Check **********************************
    stage('Git Check') {
      agent any
      steps {
        sh "git config core.whitespace -blank-at-eof"
        sh "git diff --check `git merge-base origin/master HEAD` HEAD -- . ':!*.md' ':!*.pandoc'"
        dir('scripts/jenkins') { stash(name: 'known_hosts', includes: 'known_hosts') }

        // ********* Check changesets for conditional stage execution **********
        script {
          if (currentBuild.number == 1) {
            stage_required.full = true
            return true
          }
          if (env.JOB_NAME == "ufz/ogs/master") {
            stage_required.web = true
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
                if (path.startsWith("web") && !stage_required.web) {
                  stage_required.web = true
                  echo "Doing web build."
                }
                if (path.matches("^(CMakeLists.txt|scripts|Applications|BaseLib|FileIO|GeoLib|MaterialLib|MathLib|MeshGeoToolsLib|MeshLib|NumLib|ProcessLib|SimpleTests|Tests).*")
                  && !stage_required.build) {
                  stage_required.build = true
                  echo "Doing regular build."
                }
                if (path.startsWith("Tests/Data") && !stage_required.data) {
                  stage_required.data = true
                  echo "Updating Tests/Data."
                }
              }
            }
          }
        }
      }
    }
    stage('Build') {
      parallel {
        // ************************ Docker-Conan *******************************
        stage('Docker-Conan') {
          when { expression { return stage_required.build || stage_required.full } }
          agent {
            dockerfile {
              filename 'Dockerfile.gcc.full'
              dir 'scripts/docker'
              label 'docker'
              args '-v ccache:/home/jenkins/cache/ccache -v conan-cache:/home/jenkins/cache/conan'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              // Install web dependencies
              sh("""
                cd web
                yarn --ignore-engines --non-interactive
                sudo -H pip install -r requirements.txt
                """.stripIndent())

              lock(resource: "conanCache-${env.NODE_NAME}") {
                sh 'find $CONAN_USER_HOME -name "system_reqs.txt" -exec rm {} \\;'
                configure {
                  cmakeOptions =
                    '-DOGS_USE_CONAN=ON ' +
                    '-DOGS_CONAN_BUILD=never ' +
                    '-DOGS_CPU_ARCHITECTURE=generic ' +
                    '-DDOCS_GENERATE_LOGFILE=ON ' // redirects to build/DoxygenWarnings.log
                }
              }
              build { }
              build { target="tests" }
              build { target="ctest" }
              build { target="doc" }
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
              dir('build/docs') { stash(name: 'doxygen') }
            }
            failure {
                dir('build') { deleteDir() }
            }
            success {
              script {
                publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: true,
                  keepAll: true, reportDir: 'build/docs', reportFiles: 'index.html',
                  reportName: 'Doxygen'])
                step([$class: 'WarningsPublisher', canResolveRelativePaths: false,
                  messagesPattern: """
                    .*DOT_GRAPH_MAX_NODES.
                    .*potential recursive class relation.*""",
                  parserConfigurations: [[parserName: 'Doxygen', pattern:
                  'build/DoxygenWarnings.log']], unstableTotalAll: '0'])
                archiveArtifacts 'build/*.tar.gz,build/conaninfo.txt'
                dir('build') { deleteDir() }
              }
            }
          }
        }
        // ********************* Docker-Conan-Debug ****************************
        stage('Docker-Conan-Debug') {
          when { expression { return stage_required.build || stage_required.full } }
          agent {
            dockerfile {
              filename 'Dockerfile.gcc.minimal'
              dir 'scripts/docker'
              label 'docker'
              args '-v ccache:/home/jenkins/cache/ccache -v conan-cache:/home/jenkins/cache/conan'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              lock(resource: "conanCache-${env.NODE_NAME}") {
                sh 'find $CONAN_USER_HOME -name "system_reqs.txt" -exec rm {} \\;'
                configure {
                  cmakeOptions =
                    '-DOGS_USE_CONAN=ON ' +
                    '-DOGS_CONAN_BUILD=never ' +
                    '-DOGS_CPU_ARCHITECTURE=generic '
                  config = 'Debug'
                }
              }
              build { }
              build { target = 'tests' }
            }
          }
          post {
            always {
              publishReports { }
              dir('build') { deleteDir() }
            }
          }
        }
        // ************************** envinf1 **********************************
        stage('Envinf1 (serial)') {
          when { expression { return stage_required.build || stage_required.full } }
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
          when { expression { return stage_required.build || stage_required.full } }
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
          when { expression { return stage_required.build || stage_required.full } }
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
              warnings(canResolveRelativePaths: false,
                  consoleParsers: [[parserName: 'MSBuild']],
                  excludePattern: '.*\\.conan.*',
                  messagesPattern: '.*QVTK.*')
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
        // ****************************** Mac **********************************
        stage('Mac') {
          when { expression { return stage_required.build || stage_required.full } }
          agent { label "mac"}
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DOGS_USE_CONAN=ON ' +
                  '-DOGS_CONAN_BUILD=never ' +
                  '-DOGS_CPU_ARCHITECTURE=core2 ' +
                  '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ' +
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_METIS=ON ' +
                  '-DCMAKE_OSX_DEPLOYMENT_TARGET="10.13" '
              }
              build {
                target = 'tests'
                cmd_args = '-j $(( `sysctl -n hw.ncpu` - 2 ))'
              }
              build {
                target = 'ctest'
                cmd_args = '-j $(( `sysctl -n hw.ncpu` - 2 ))'
              }
              build {
                target = 'package'
                cmd_args = '-j $(( `sysctl -n hw.ncpu` - 2 ))'
              }
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
                archiveArtifacts 'build/*.tar.gz,build/*.dmg,build/conaninfo.txt'
                dir('build') { deleteDir() }
            }
          }
        }
        // **************************** Web ************************************
        stage('Web') {
          when { expression { return stage_required.web || stage_required.full } }
          agent {
            dockerfile {
              filename 'Dockerfile.gcc.full'
              dir 'scripts/docker'
              label 'docker'
              additionalBuildArgs '--pull'
            }
          }
          environment {
            CONTENTFUL_ACCESS_TOKEN = credentials('CONTENTFUL_ACCESS_TOKEN')
            CONTENTFUL_OGS_SPACE_ID = credentials('CONTENTFUL_OGS_SPACE_ID')
            ALGOLIA_WRITE_KEY = credentials('ALGOLIA_WRITE_KEY')
          }
          steps {
            dir ('web') {
              sh "yarn --ignore-engines --ignore-optional --non-interactive"
              sh "pandoc-citeproc --bib2json ../Documentation/bibliography.bib > data/bibliography.json"
              sh "node_modules/.bin/webpack -p"
              script {
                if (env.JOB_NAME == 'ufz/ogs/master') {
                  sh "hugo --ignoreCache --baseURL https://benchmarks.opengeosys.org"
                  sh ("node_modules/.bin/hugo-algolia --toml -s")
                } else {
                  sh ("hugo --ignoreCache --baseURL " + env.JOB_URL + "Web/")
                }
              }
            }
          }
          post {
            success {
              dir('web/public') { stash(name: 'web') }
              script {
                publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: true,
                  keepAll: true, reportDir: 'web/public', reportFiles: 'index.html',
                  reportName: 'Web'])
              }
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
        // ************************* Analyzers *********************************
        stage('Analyzers') {
          when { expression { return stage_required.build || stage_required.full } }
          agent {
            dockerfile {
              filename 'Dockerfile.clang.full'
              dir 'scripts/docker'
              label 'docker'
              args '-v ccache:/home/jenkins/cache/ccache -v conan-cache:/home/jenkins/cache/conan'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              lock(resource: "conanCache-${env.NODE_NAME}") {
                sh 'find $CONAN_USER_HOME -name "system_reqs.txt" -exec rm {} \\;'
                configure {
                  cmakeOptions =
                    '-DOGS_USE_CONAN=ON ' +
                    '-DOGS_CONAN_BUILD=never ' +
                    '"-DCMAKE_CXX_INCLUDE_WHAT_YOU_USE=include-what-you-use;-Xiwyu;--mapping_file=../scripts/jenkins/iwyu-mappings.imp" ' +
                    '-DCMAKE_LINK_WHAT_YOU_USE=ON ' +
                    '"-DCMAKE_CXX_CPPCHECK=cppcheck;--std=c++11;--language=c++;--suppress=syntaxError;--suppress=preprocessorErrorDirective:*/ThirdParty/*;--suppress=preprocessorErrorDirective:*conan*/package/*" ' +
                    '-DCMAKE_CXX_CLANG_TIDY=clang-tidy-3.9 '
                  config = 'Release'
                }
              }
              build { target = 'check-header' }
              build { }
            }
          }
          post { always { dir('build') { deleteDir() } } }
        }
        // ************************* Deploy Web ********************************
        stage('Deploy Web') {
          when { expression { return stage_required.web || stage_required.full } }
          agent any
          steps {
            dir('web') { unstash 'web' }
            unstash 'known_hosts'
            script {
              sshagent(credentials: ['www-data_jenkins']) {
                sh 'rsync -a --delete --stats -e "ssh -o UserKnownHostsFile=' +
                   'known_hosts" web/. ' +
                   'www-data@jenkins.opengeosys.org:/var/www/dev.opengeosys.org'
              }
            }
          }
        }
        stage('Deploy Doxygen') {
          when { expression { return stage_required.build || stage_required.full } }
          agent any
          steps {
            dir('doxygen') { unstash 'doxygen' }
            unstash 'known_hosts'
            script {
              sshagent(credentials: ['www-data_jenkins']) {
                sh 'rsync -a --delete --stats -e "ssh -o UserKnownHostsFile=' +
                   'known_hosts" doxygen/. ' +
                   'www-data@jenkins.opengeosys.org:/var/www/doxygen.opengeosys.org'
              }
            }
          }
        }
        // *********************** Deploy envinf1 ******************************
        stage('Deploy envinf1') {
          when { expression { return stage_required.build || stage_required.full } }
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
          when { expression { return stage_required.build || stage_required.full } }
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
          when { expression { return stage_required.build || stage_required.full } }
          agent {
            dockerfile {
              filename 'Dockerfile.clang.minimal'
              dir 'scripts/docker'
              label 'docker'
              args '-v ccache:/home/jenkins/cache/ccache -v conan-cache:/home/jenkins/cache/conan'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              sh 'find $CONAN_USER_HOME -name "system_reqs.txt" -exec rm {} \\;'
              configure {
                cmakeOptions =
                  '-DOGS_USE_CONAN=ON ' +
                  '-DOGS_ADDRESS_SANITIZER=ON ' +
                  '-DOGS_UNDEFINED_BEHAVIOR_SANITIZER=ON ' +
                  '-DOGS_BUILD_UTILS=ON '
              }
              try {
                build { cmd = 'UBSAN_OPTIONS=print_stacktrace=1 ninja tests' }
              }
              catch(err) { echo "Clang sanitizer for unit tests failed!" }

              try {
                build { cmd = 'UBSAN_OPTIONS=print_stacktrace=1 ninja ctest' }
              }
              catch(err) { echo "Clang sanitizer for end-to-end tests failed!" }
            }
          }
          post {
            always {
              dir('build') { deleteDir() }
              warnings(canResolveRelativePaths: false,
                  consoleParsers: [[parserName: 'Clang (LLVM based)']])
            }
          }
        }
        // ********************* Update ufz/ogs-data ***************************
        stage('Update ogs-data') {
          when { expression { return stage_required.data } }
          agent any
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
      // post {
      //   always {
      //     step([$class: 'AnalysisPublisher', unstableNewAll: '1'])
      //   }
      // }
    } // end stage master
  }
}
