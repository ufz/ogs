#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.12') _

def stage_required = [build: false, data: false, full: false, docker: false]

pipeline {
  agent none
  options {
    ansiColor('xterm')
    timestamps()
    buildDiscarder(logRotator(numToKeepStr: '30', artifactNumToKeepStr: '10'))
    timeout(time: 3, unit: 'HOURS')
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
          def causes = currentBuild.rawBuild.getCauses()
          for(cause in causes) {
            if (cause.class.toString().contains("UserIdCause")) {
              echo "Doing full build because job was started by user."
              stage_required.full = true
              return true
            }
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
                if (path.startsWith("Tests/Data") && !stage_required.data) {
                  stage_required.data = true
                  echo "Updating Tests/Data."
                }
                if (path.startsWith("scripts/docker") && !stage_required.docker) {
                  stage_required.docker = true
                  echo "Doing Docker images build."
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
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.gcc.full'
              dir 'scripts/docker'
              label 'docker'
              args '-v /datadrive/cache:/home/jenkins/cache'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              sh 'conan user'
              configure {
                cmakeOptions =
                  '-DOGS_CPU_ARCHITECTURE=generic ' +
                  '-DDOCS_GENERATE_LOGFILE=ON ' + // redirects to build/DoxygenWarnings.log
                  '-DOGS_USE_PYTHON=ON '
              }
              build { }
              build { target="tests" }
              build { target="ctest" }
              build { target="doc" }
              // TODO: .*DOT_GRAPH_MAX_NODES.
              //       .*potential recursive class relation.*
              recordIssues tools: [[pattern: 'build/DoxygenWarnings.log',
                tool: [$class: 'Doxygen']]],
                unstableTotalAll: 23
              dir('build/docs') { stash(name: 'doxygen') }
              publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: true,
                  keepAll: true, reportDir: 'build/docs', reportFiles: 'index.html',
                  reportName: 'Doxygen'])
              configure {
                cmakeOptions =
                  '-DOGS_CPU_ARCHITECTURE=generic ' +
                  '-DOGS_USE_PYTHON=ON ' +
                  '-DOGS_USE_PCH=OFF ' +     // see #1992
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_TESTS=OFF '
              }
              build { log="build.log" }
              build { target="doc" }
            }
          }
          post {
            always {
              publishReports { }
              recordIssues enabledForFailure: true, filters: [
                excludeFile('.*qrc_icons\\.cpp.*'), excludeFile('.*QVTKWidget.*')],
                tools: [[pattern: 'build/build.log', name: 'GCC',
                  tool: [$class: 'GnuMakeGcc']]],
                unstableTotalAll: 19
            }
            success {
              dir('build/docs') { stash(name: 'doxygen') }
              archiveArtifacts 'build/*.tar.gz,build/conaninfo.txt'
            }
          }
        }
        // ********************* Docker-Conan-Debug ****************************
        stage('Docker-Conan-Debug') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.gcc.minimal'
              dir 'scripts/docker'
              label 'docker'
              args '-v /datadrive/cache:/home/jenkins/cache'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              sh 'conan user'
              configure {
                cmakeOptions =
                  '-DOGS_CPU_ARCHITECTURE=generic '
                config = 'Debug'
              }
              build { }
              build { target = 'tests' }
            }
          }
          post {
            always { publishReports { } }
          }
        }
        // ************************** envinf1 **********************************
        stage('Envinf1 (serial)') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent { label "envinf1"}
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DOGS_BUILD_UTILS=ON ' +
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
            always { publishReports { } }
          }
        }
        stage('Envinf1 (parallel)') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent { label "envinf1"}
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DOGS_BUILD_UTILS=ON ' +
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
            always { publishReports { } }
          }
        }
        // ************************** Windows **********************************
        stage('Win') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent {label 'win && conan' }
          environment {
            MSVC_NUMBER = '15'
            MSVC_VERSION = '2017'
          }
          steps {
            script {
              // CLI
              bat 'conan remove --locks'
              configure {
                cmakeOptions =
                  '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ' +
                  '-DOGS_USE_PYTHON=ON '
              }
              build { }
              build { target="tests" }
              build { target="ctest" }
              // GUI
              configure {
                cmakeOptions =
                  '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ' +
                  '-DOGS_USE_PYTHON=ON ' +
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_TESTS=OFF ' +
                  '-DOGS_BUILD_SWMM=ON '
              }
              build { log="build.log" }
            }
          }
          post {
            always {
              publishReports { }
              recordIssues enabledForFailure: true, filters: [
                excludeFile('.*\\.conan.*'), excludeFile('.*ThirdParty.*'),
                excludeFile('.*thread.hpp')],
                tools: [[pattern: 'build/build.log', name: 'MSVC',
                  tool: [$class: 'MsBuild']]],
                unstableTotalAll: 4
            }
            success {
              archiveArtifacts 'build/*.zip,build/conaninfo.txt'
            }
          }
        }
        // ****************************** Mac **********************************
        stage('Mac') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent { label "mac"}
          steps {
            script {
              configure {
                cmakeOptions =
                  '-DOGS_CPU_ARCHITECTURE=core2 ' +
                  '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ' +
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DCMAKE_OSX_DEPLOYMENT_TARGET="10.13" '
              }
              build {
                log = "build.log"
                cmd_args = '-j $(( `sysctl -n hw.ncpu` - 2 ))'
              }
              build {
                target = 'tests'
                cmd_args = '-j $(( `sysctl -n hw.ncpu` - 2 ))'
              }
              build {
                target = 'ctest'
                cmd_args = '-j $(( `sysctl -n hw.ncpu` - 2 ))'
              }
            }
          }
          post {
            always {
              publishReports { }
              recordIssues enabledForFailure: true, filters: [
                excludeFile('.*qrc_icons\\.cpp.*'), excludeFile('.*QVTKWidget.*')],
                tools: [[pattern: 'build/build.log', name: 'Clang (macOS)',
                   id: 'clang-mac', tool: [$class: 'Clang']]],
                unstableTotalAll: 3
            }
            success {
              archiveArtifacts 'build/*.tar.gz,build/*.dmg,build/conaninfo.txt'
            }
          }
        }
      } // end parallel
    } // end stage Build
    stage('Master') {
      when { environment name: 'JOB_NAME', value: 'ufz/ogs/master' }
      parallel {
        // ********************* Push Docker Images ****************************
        stage('Push Docker Images') {
          when {
            beforeAgent true
            expression { return stage_required.docker || stage_required.full }
          }
          agent { label 'docker' }
          steps {
            script {
              dir('scripts/docker') {
                def gccImage = docker.build("ogs6/gcc:latest", "-f Dockerfile.gcc.full .")
                def clangImage = docker.build("ogs6/clang:latest", "-f Dockerfile.clang.full .")
                docker.withRegistry('https://registry.hub.docker.com', 'docker-hub-credentials') {
                  gccImage.push()
                  clangImage.push()
                }
              }
            }
          }
        }
        // ************************* Analyzers *********************************
        stage('Analyzers') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.clang.full'
              dir 'scripts/docker'
              label 'docker'
              args '-v /datadrive/cache:/home/jenkins/cache'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              lock(resource: "conanCache-${env.NODE_NAME}") {
                sh 'conan user'
                sh 'find $CONAN_USER_HOME -name "system_reqs.txt" -exec rm {} \\;'
                configure {
                  cmakeOptions =
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
        }
        // *********************** Deploy Doxygen ******************************
        stage('Deploy Doxygen') {
          when { expression { return stage_required.build || stage_required.full } }
          agent { label "master" }
          steps {
            dir('doxygen') { unstash 'doxygen' }
            unstash 'known_hosts'
            script {
              sshagent(credentials: ['www-data_jenkins']) {
                sh 'rsync -a --delete --stats -e "ssh -o UserKnownHostsFile=' +
                   'known_hosts" doxygen/. ' +
                   'www-data@jenkins:/var/www/doxygen.opengeosys.org'
              }
            }
          }
        }
        // *********************** Deploy envinf1 ******************************
        stage('Deploy envinf1') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent { label "envinf1"}
          steps {
            script {
              sh 'rm -rf /global/apps/ogs/head/standard'
              configure {
                cmakeOptions =
                  '-DOGS_BUILD_UTILS=ON ' +
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
        }
        // ******************** Deploy envinf1 PETSc ***************************
        stage('Deploy envinf1 PETSc') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent { label "envinf1"}
          steps {
            script {
              sh 'rm -rf /global/apps/ogs/head/petsc'
              configure {
                cmakeOptions =
                  '-DOGS_USE_PETSC=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
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
        }
        // ************************** Sanitizer ********************************
        stage('Sanitizer') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.clang.minimal'
              dir 'scripts/docker'
              label 'docker'
              args '-v /datadrive/cache:/home/jenkins/cache'
              additionalBuildArgs '--pull'
            }
          }
          environment {
            UBSAN_OPTIONS = 'print_stacktrace=1'
            LSAN_OPTIONS = "suppressions=$WORKSPACE/scripts/test/leak_sanitizer.suppressions"
          }
          steps {
            script {
              sh 'conan user'
              sh 'find $CONAN_USER_HOME -name "system_reqs.txt" -exec rm {} \\;'
              configure {
                cmakeOptions =
                  '-DOGS_ADDRESS_SANITIZER=ON ' +
                  '-DOGS_UNDEFINED_BEHAVIOR_SANITIZER=ON ' +
                  '-DOGS_BUILD_UTILS=ON '
              }
              try {
                build { target = 'tests' }
              }
              catch(err) { echo "Clang sanitizer for unit tests failed!" }

              try {
                build { target = 'ctest' }
              }
              catch(err) { echo "Clang sanitizer for end-to-end tests failed!" }
            }
          }
          post {
            always {
              recordIssues enabledForFailure : true,
                filters: [includeCategory('clang-analyzer.*')],
                tools: [[name:'Clang (StaticAnalyzer)', tool:[$class:'Clang']]]
            }
          }
        }
        // ********************* Update ufz/ogs-data ***************************
        stage('Update ogs-data') {
          when {
            beforeAgent true
            expression { return stage_required.data }
          }
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
    } // end stage master
  }
}
