#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.19') _

def stage_required = [build: false, data: false, full: false, docker: false]
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
    booleanParam(name: 'envinf1_serial', defaultValue: true)
    booleanParam(name: 'envinf1_parallel', defaultValue: true)
    booleanParam(name: 'win', defaultValue: true)
    booleanParam(name: 'mac', defaultValue: true)
    booleanParam(name: 'clang_analyzer', defaultValue: true)
  }
  stages {
     // *************************** Git Check **********************************
    stage('Git Check') {
      agent { label "master"}
      steps {
        sh "git config core.whitespace -blank-at-eof"
        sh "git diff --check `git merge-base origin/master HEAD` HEAD -- . ':!*.md' ':!*.pandoc' ':!*.asc'"
        dir('scripts/jenkins') { stash(name: 'known_hosts', includes: 'known_hosts') }
        ciSkip action: 'check' // Check for [ci skip] commit message.

        // ********* Check changesets for conditional stage execution **********
        script {
          if (env.JOB_NAME == 'ufz/ogs/master') {
            build_shared = 'OFF'
          }
          if (currentBuild.number == 1 || buildingTag()) {
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
          }
          steps {
            script {
              sh 'git submodule sync'
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=${build_shared} " +
                  '-DOGS_CPU_ARCHITECTURE=generic ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_USE_CVODE=ON '
              }
              build {
                target="package"
                log="build1.log"
              }
              configure { // CLI + Python
                cmakeOptions = '-DOGS_USE_PYTHON=ON '
                keepDir = true
              }
              build {
                target="package"
                log="build2.log"
              }
              build { target="tests" }
              build { target="ctest" }
              build { target="doc" }
            }
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
                failedTotalAll: 1
            }
            success {
              publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: true,
                  keepAll: true, reportDir: 'build/docs', reportFiles: 'index.html',
                  reportName: 'Doxygen'])
              archiveArtifacts 'build/*.tar.gz,build/conaninfo.txt'
              dir('build/docs') { stash(name: 'doxygen') }
            }
          }
        }
        // *********************** Docker-Conan-GUI *****************************
        stage('Docker-Conan-GUI') {
          when {
            beforeAgent true
            expression { return params.docker_conan_gui && (stage_required.build || stage_required.full) }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.gcc.gui'
              dir 'scripts/docker'
              // Singularity1 has on old kernel (3.10) which is not compatible with Qt > 5.10 (req. 3.15)
              label 'docker && !singularity1'
              args '-v /home/jenkins/cache/ccache:/opt/ccache -v /home/jenkins/cache/conan/.conan:/opt/conan/.conan'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              sh 'git submodule sync'
              sh 'conan remove --system-reqs Qt/5.11.2@bilke/stable'
              sh 'conan remove --system-reqs VTK/8.1.1@bilke/stable'
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=${build_shared} " +
                  '-DOGS_CPU_ARCHITECTURE=generic ' +
                  '-DOGS_USE_PCH=OFF ' +     // see #1992
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_TESTS=OFF '
              }
              build {
                target="package"
                log="build1.log"
              }
              configure { // CLI + GUI + Python
                cmakeOptions = '-DOGS_USE_PYTHON=ON '
                keepDir = true
              }
              build {
                target="package"
                log="build2.log"
              }
              build { target="cppcheck" }
            }
          }
          post {
            always {
              recordIssues enabledForFailure: true, filters: [
                excludeFile('.*qrc_icons\\.cpp.*'),
                excludeMessage('.*tmpnam.*')],
                tools: [gcc4(name: 'GCC-GUI', id: 'gcc4-gui',
                             pattern: 'build/build*.log')],
                unstableTotalAll: 1
              recordIssues enabledForFailure: true,
                tools: [cppCheck(pattern: 'build/cppcheck.log')]
            }
            success { archiveArtifacts 'build/*.tar.gz,build/conaninfo.txt' }
          }
        }
        // ********************* Docker-Conan-Debug ****************************
        stage('Docker-Conan-Debug') {
          when {
            beforeAgent true
            expression { return params.docker_conan_debug && (stage_required.build || stage_required.full) }
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
          steps {
            script {
              sh 'git submodule sync'
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=${build_shared} " +
                  '-DOGS_CPU_ARCHITECTURE=generic ' +
                  '-DOGS_COVERAGE=ON '
                config = 'Debug'
              }
              build { }
              build { target = 'testrunner_coverage_cobertura' }
              dir('build') {
                sh "curl -s https://codecov.io/bash | bash -s - -f *_cobertura.xml"
              }
            }
          }
          post {
            always {
              xunit([GoogleTest(pattern: 'build/Tests/testrunner.xml')])
            }
          }
        }
        // ************************** envinf1 **********************************
        stage('Envinf1 (serial)') {
          when {
            beforeAgent true
            expression { return params.envinf1_serial && (stage_required.build || stage_required.full) }
          }
          agent { label "envinf1"}
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
                  '-DOGS_USE_CONAN=OFF '
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
              xunit([
                CTest(pattern: 'build/Testing/**/*.xml'),
                GoogleTest(pattern: 'build/Tests/testrunner.xml')
              ])
            }
          }
        }
        stage('Envinf1 (parallel)') {
          when {
            beforeAgent true
            expression { return params.envinf1_parallel && (stage_required.build || stage_required.full) }
          }
          agent { label "envinf1"}
          environment {
            OMP_NUM_THREADS = '1'
          }
          steps {
            script {
              configure {
                sh 'git submodule sync'
                cmakeOptions =
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DBUILD_SHARED_LIBS=ON ' +
                  '-DOGS_USE_PETSC=ON ' +
                  '-DOGS_USE_CONAN=OFF '
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
              xunit([
                CTest(pattern: 'build/Testing/**/*.xml'),
                GoogleTest(pattern: 'build/Tests/testrunner.xml')
              ])
            }
          }
        }
        // ************************** Windows **********************************
        stage('Win') {
          when {
            beforeAgent true
            expression { return params.win && (stage_required.build || stage_required.full) }
          }
          agent {label 'win && conan' }
          environment {
            MSVC_NUMBER = '15'
            MSVC_VERSION = '2017'
            OMP_NUM_THREADS = '1'
          }
          steps {
            script {
              def num_threads = env.NUM_THREADS
              bat 'git submodule sync'
              bat 'conan remove --locks'
              configure { // CLI + GUI
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=OFF " +
                  '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ' +
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_SWMM=ON '
              }
              build {
                target="package"
                log="build1.log"
                cmd_args="-l ${num_threads}"
              }
              configure { // CLI + GUI + Python
                cmakeOptions = '-DOGS_USE_PYTHON=ON '
                keepDir = true
              }
              build {
                target="package"
                log="build2.log"
              }
              build { target="tests" }
              build { target="ctest" }
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
        // ****************************** Mac **********************************
        stage('Mac') {
          when {
            beforeAgent true
            expression { return params.mac && (stage_required.build || stage_required.full) }
          }
          agent { label "mac"}
          environment {
            OMP_NUM_THREADS = '1'
          }
          steps {
            script {
              sh 'git submodule sync'
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=${build_shared} " +
                  '-DOGS_CPU_ARCHITECTURE=core2 ' +
                  '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ' +
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DCMAKE_OSX_DEPLOYMENT_TARGET="10.13" '
              }
              build {
                target="package"
                log = "build.log"
              }
              build { target = 'tests' }
              build { target = 'ctest' }
            }
          }
          post {
            always {
              xunit([
                CTest(pattern: 'build/Testing/**/*.xml'),
                GoogleTest(pattern: 'build/Tests/testrunner.xml')
              ])
              recordIssues enabledForFailure: true, filters: [
                excludeFile('.*qrc_icons\\.cpp.*'), excludeMessage('.*QVTKWidget.*'),
                excludeMessage('.*tmpnam.*')],
                tools: [clang(name: 'Clang (macOS)', pattern: 'build/build.log',
                  id: 'clang-mac')], unstableTotalAll: 1
            }
            success {
              archiveArtifacts 'build/*.tar.gz,build/*.dmg,build/conaninfo.txt'
            }
          }
        }
        // ************************* Clang-Analyzer *********************************
        stage('Clang-Analyzer') {
          when {
            beforeAgent true
            expression { return params.clang_analyzer && (stage_required.build || stage_required.full) }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.clang.full'
              dir 'scripts/docker'
              label 'docker'
              args '-v /home/jenkins/cache/ccache:/opt/ccache -v /home/jenkins/cache/conan/.conan:/opt/conan/.conan'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              sh 'git submodule sync'
              sh 'find $CONAN_USER_HOME -name "system_reqs.txt" -exec rm {} \\;'
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=${build_shared} " +
                  '-DBUILD_TESTING=OFF ' +
                  '-DCMAKE_CXX_CLANG_TIDY=clang-tidy-7 '
              }
              build { log = 'build.log' }
            }
          }
          post {
            always {
              recordIssues enabledForFailure: true, filters: [
                excludeFile('.*\\.conan.*')],
                tools: [clangTidy(name: 'Clang-Tidy', pattern: 'build/build.log')],
                qualityGates: [[threshold: 513, type: 'TOTAL', unstable: true]]
            }
          }
        }
      } // end parallel
    } // end stage Build
    stage('Master') {
      when { environment name: 'JOB_NAME', value: 'ufz/ogs/master' }
      parallel {
        // ************************* Tests-Large *******************************
        stage('Tests-Large') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.gcc.full'
              dir 'scripts/docker'
              label 'envinf11w || envinf56'
              args '-v /home/jenkins/cache/ccache:/opt/ccache -v /home/jenkins/cache/conan/.conan:/opt/conan/.conan'
              additionalBuildArgs '--pull'
            }
          }
          environment {
            OMP_NUM_THREADS = '1'
          }
          steps {
            script {
              configure { }
              build { target = 'ctest-large' }
            }
          }
          post {
            always {
              xunit([CTest(pattern: 'build/Testing/**/*.xml')])
              dir('build') { deleteDir() }
            }
          }
        }
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
                def gccGuiImage = docker.build("ogs6/gcc:gui", "-f Dockerfile.gcc.gui .")
                def clangImage = docker.build("ogs6/clang:latest", "-f Dockerfile.clang.full .")
                docker.withRegistry('https://registry.hub.docker.com', 'docker-hub-credentials') {
                  gccImage.push()
                  gccGuiImage.push()
                  clangImage.push()
                }
              }
            }
          }
        }
        // ************************* Check headers *********************************
        stage('Check headers') {
          when {
            beforeAgent true
            expression { return stage_required.build || stage_required.full }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.clang.full'
              dir 'scripts/docker'
              label 'docker'
              args '-v /home/jenkins/cache/ccache:/opt/ccache -v /home/jenkins/cache/conan/.conan:/opt/conan/.conan'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              sh 'git submodule sync'
              sh 'find $CONAN_USER_HOME -name "system_reqs.txt" -exec rm {} \\;'
              try {
                configure {
                  cmakeOptions = '-DOGS_CHECK_HEADER_COMPILATION=ON'
                  dir = 'build-check-header'
                }
              }
              catch(err) {
                echo "check-header failed!"
                sh 'cat build-check-header/CMakeFiles/CMakeError.log'
              }
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
              args '-v /home/jenkins/cache/ccache:/opt/ccache -v /home/jenkins/cache/conan/.conan:/opt/conan/.conan'
              additionalBuildArgs '--pull'
            }
          }
          environment {
            UBSAN_OPTIONS = 'print_stacktrace=1'
            LSAN_OPTIONS = "suppressions=$WORKSPACE/scripts/test/leak_sanitizer.suppressions"
          }
          steps {
            script {
              sh 'git submodule sync'
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
          // Currently disabled because of Java out ouf heap space errors
          // post {
            // always {
              // recordIssues enabledForFailure : true,
                // filters: [includeCategory('clang-analyzer.*')],
                // tools: [clang(name: 'Clang (StaticAnalyzer)')]
            // }
          // }
        }
        // ********************* Update ufz/ogs-data ***************************
        stage('Update ogs-data') {
          when {
            beforeAgent true
            expression { return stage_required.data }
          }
          agent { label "master"}
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
          agent { label "master"}
          when { buildingTag() }
          steps {
            keepBuild()
            script { currentBuild.displayName = tag }
          }
        }
      } // end parallel
    } // end stage master
  }
}
