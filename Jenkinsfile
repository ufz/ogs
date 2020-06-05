#!/usr/bin/env groovy
@Library('jenkins-pipeline@1.0.23') _

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
    booleanParam(name: 'mac_gui', defaultValue: true)
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
              // Workaround some MGIS CMake logic flaws
              configure {
                keepDir = true
              }
              build {
                target="package"
                log="build1.log"
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
                excludeFile('.*\\.conan.*/TFEL/.*'),
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
                unstableTotalAll: 2
            }
            success {
              publishHTML(target: [allowMissing: false, alwaysLinkToLastBuild: true,
                  keepAll: true, reportDir: 'build/docs', reportFiles: 'index.html',
                  reportName: 'Doxygen'])
              archiveArtifacts 'build/*.tar.gz,build/conaninfo.txt'
              script {
                if (env.JOB_NAME == 'ufz/ogs/master') {
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
              label 'docker'
              args '-v /home/jenkins/cache/ccache:/opt/ccache -v /home/jenkins/cache/conan/.conan:/opt/conan/.conan'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              sh 'git submodule sync && git submodule update'
              sh "conan remove --system-reqs '*'"
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=${build_shared} " +
                  '-DOGS_CPU_ARCHITECTURE=generic ' +
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_CONAN_BUILD=missing ' +
                  '-DOGS_BUILD_TESTS=OFF ' +
                  '-DOGS_USE_NETCDF=ON '
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
                excludeFile('ncGroup\\.h'),
                excludeMessage('.*tmpnam.*')],
                tools: [gcc4(name: 'GCC-GUI', id: 'gcc4-gui',
                             pattern: 'build/build*.log')],
                unstableTotalAll: 1
              recordIssues enabledForFailure: true,
                tools: [cppCheck(pattern: 'build/cppcheck.log')],
                unstableTotalAll: 500
            }
            success { archiveArtifacts 'build/*.tar.gz,build/conaninfo.txt,build/cppcheck.log' }
          }
        }
        // ********************* Docker-Conan-Debug ****************************
        stage('Docker-Conan-Debug') {
          when {
            beforeAgent true
            // expression { return params.docker_conan_debug && (stage_required.build || stage_required.full) }
            expression { return false }
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
              sh 'git submodule sync && git submodule update'
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=${build_shared} " +
                  '-DOGS_CONAN_BUILD=missing ' +
                  '-DOGS_CONAN_BUILD_TYPE=Release ' +
                  '-DOGS_CPU_ARCHITECTURE=generic '
                config = 'Debug'
              }
              build { }
              build { target = 'tests' }
            }
          }
          post {
            always {
              xunit([GoogleTest(pattern: 'build/Tests/testrunner.xml')])
            }
          }
        }
        // **************************** eve ************************************
        stage('Frontend2 (serial)') {
          when {
            beforeAgent true
            expression { return params.eve_serial && (stage_required.build || stage_required.full) }
          }
          agent { label "frontend2"}
          environment {
            OMP_NUM_THREADS = '1'
            SOURCE_DIR = "${env.WORKSPACE}"
            BUILD_DIR = "/tmp/${env.BUILD_TAG}"
          }
          steps {
            script {
              sh 'git submodule sync && git submodule update'
              configure {
                cmakeOptions =
                  '-DOGS_USE_CONAN=OFF ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DBUILD_SHARED_LIBS=ON ' +
                  '-DOGS_CPU_ARCHITECTURE=sandybridge ' +
                  '-DCMAKE_INSTALL_PREFIX=/global/apps/ogs/head/standard ' +
                  '-DOGS_MODULEFILE=/global/apps/modulefiles/ogs/head/standard '
                env = 'eve/cli.sh'
              }
              build {
                env = 'eve/cli.sh'
                cmd = 'nice -n 15 cmake --build . --config Release'
                cmd_args = '-j 8'
              }
              build {
                env = 'eve/cli.sh'
                cmd = 'nice -n 15 cmake --build . --config Release --target tests'
              }
              build {
                env = 'eve/cli.sh'
                cmd = 'nice -n 15 cmake --build . --config Release --target ctest'
              }
            }
          }
          post {
            success {
              script {
                if (env.JOB_NAME == 'ufz/ogs/master') {
                  sh 'rm -rf /global/apps/ogs/head/standard'
                  build {
                    env = 'eve/cli.sh'
                    cmd = 'nice -n 15 cmake --build . --config Release --target install'
                  }
                }
              }
            }
            always {
              sh "mkdir _out && cp -r ${env.BUILD_DIR}/Testing _out/ && cp -r ${env.BUILD_DIR}/Tests/testrunner.xml _out/"
              xunit([
                CTest(pattern: "_out/Testing/**/*.xml"),
                GoogleTest(pattern: "_out/testrunner.xml")
              ])
            }
            cleanup {
              dir("${env.BUILD_DIR}") { deleteDir() }
              dir('_out') { deleteDir() }
            }
          }
        }
        stage('Frontend2 (parallel)') {
          when {
            beforeAgent true
            expression { return params.eve_parallel && (stage_required.build || stage_required.full) }
          }
          agent { label "frontend2"}
          environment {
            OMP_NUM_THREADS = '1'
            SOURCE_DIR = "${env.WORKSPACE}"
            BUILD_DIR = "/tmp/${env.BUILD_TAG}-petsc"
          }
          steps {
            script {
              sh 'git submodule sync && git submodule update'
              configure {
                cmakeOptions =
                  '-DOGS_USE_CONAN=OFF ' +
                  '-DOGS_USE_PETSC=ON ' +
                  '-DBUILD_SHARED_LIBS=ON ' +
                  '-DOGS_CPU_ARCHITECTURE=sandybridge ' +
                  '-DCMAKE_INSTALL_PREFIX=/global/apps/ogs/head/petsc ' +
                  '-DOGS_MODULEFILE=/global/apps/modulefiles/ogs/head/petsc '
                env = 'eve/petsc.sh'
              }
              build {
                env = 'eve/petsc.sh'
                cmd_args = '-j 8'
                cmd = 'nice -n 15 cmake --build . --config Release'
              }
              build {
                env = 'eve/petsc.sh'
                target = 'tests'
                cmd = 'nice -n 15 cmake --build . --config Release --target tests'
              }
              build {
                env = 'eve/petsc.sh'
                target = 'ctest'
                cmd = 'nice -n 15 cmake --build . --config Release --target ctest'
              }
            }
          }
          post {
            success {
              script {
                if (env.JOB_NAME == 'ufz/ogs/master') {
                  sh 'rm -rf /global/apps/ogs/head/petsc'
                  build {
                    env = 'eve/petsc.sh'
                    cmd = 'nice -n 15 cmake --build . --config Release --target install'
                  }
                }
              }
            }
            always {
              sh "mkdir _out && cp -r ${env.BUILD_DIR}/Testing _out/ && cp -r ${env.BUILD_DIR}/Tests/testrunner.xml _out/"
              xunit([
                CTest(pattern: "_out/Testing/**/*.xml"),
                GoogleTest(pattern: "_out/testrunner.xml")
              ])
            }
            cleanup {
              dir("${env.BUILD_DIR}") { deleteDir() }
              dir('_out') { deleteDir() }
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
            MSVC_NUMBER = '16'
            MSVC_VERSION = '2019'
            OMP_NUM_THREADS = '1'
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
                excludeFile('.*\\.conan.*'),
                excludeFile('.*thread.hpp')],
                tools: [msBuild(name: 'MSVC', pattern: 'build/build*.log')],
                qualityGates: [[threshold: 7, type: 'TOTAL', unstable: true]]
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
              sh 'git submodule sync && git submodule update'
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=OFF " +
                  '-DOGS_CPU_ARCHITECTURE=core2 ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_CONAN_BUILD=missing ' +
                  '-DCMAKE_OSX_DEPLOYMENT_TARGET="10.15" '
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
                  id: 'clang-mac')], unstableTotalAll: 3
            }
            success {
              archiveArtifacts 'build/*.tar.gz,build/*.dmg,build/conaninfo.txt'
            }
          }
        }
        // **************************** Mac-Gui ********************************
        stage('Mac-Gui') {
          when {
            beforeAgent true
            expression { return params.mac_gui && (stage_required.build || stage_required.full) }
          }
          agent { label "mac"}
          environment {
            OMP_NUM_THREADS = '1'
          }
          steps {
            script {
              sh 'git submodule sync && git submodule update'
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=${build_shared} " +
                  '-DOGS_CPU_ARCHITECTURE=core2 ' +
                  '-DOGS_BUILD_UTILS=ON ' +
                  '-DOGS_BUILD_GUI=ON ' +
                  '-DOGS_USE_CONAN=OFF ' +
                  '-DCMAKE_OSX_DEPLOYMENT_TARGET="10.15" ' +
                  '-DOGS_USE_NETCDF=ON '
              }
              build {
                target="package"
                log = "build.log"
              }
              build { target = 'tests' }
            }
          }
          post {
            always {
              xunit([
                GoogleTest(pattern: 'build/Tests/testrunner.xml')
              ])
              recordIssues enabledForFailure: true, filters: [
                excludeFile('.*qrc_icons\\.cpp.*'), excludeMessage('.*QVTKWidget.*'),
                excludeMessage('.*tmpnam.*')],
                tools: [clang(name: 'Clang (macOS, GUI)', pattern: 'build/build.log',
                  id: 'clang-mac-gui')], unstableTotalAll: 1
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
              sh 'git submodule sync && git submodule update'
              sh 'find $CONAN_USER_HOME -name "system_reqs.txt" -exec rm {} \\;'
              configure {
                cmakeOptions =
                  "-DBUILD_SHARED_LIBS=${build_shared} " +
                  '-DBUILD_TESTING=OFF ' +
                  '-DCMAKE_CXX_CLANG_TIDY=clang-tidy-9 '
              }
              build { log = 'build.log' }
            }
          }
          post {
            always {
              recordIssues enabledForFailure: true, filters: [
                excludeFile('.*\\.conan.*')],
                tools: [clangTidy(name: 'Clang-Tidy', pattern: 'build/build.log')],
                qualityGates: [[threshold: 165, type: 'TOTAL', unstable: true]]
            }
          }
        }
        // ************************* Check headers *********************************
        stage('Check headers') {
          when {
            beforeAgent true
            allOf {
              expression { return params.master_jobs && (stage_required.build || stage_required.full) }
              environment name: 'JOB_NAME', value: 'ufz/ogs/master'
            }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.clang.gui'
              dir 'scripts/docker'
              label 'docker'
              args '-v /home/jenkins/cache/ccache:/opt/ccache -v /home/jenkins/cache/conan/.conan:/opt/conan/.conan'
              additionalBuildArgs '--pull'
            }
          }
          steps {
            script {
              sh 'git submodule sync && git submodule update'
              sh 'find $CONAN_USER_HOME -name "system_reqs.txt" -exec rm {} \\;'
              try {
                configure {
                  cmakeOptions =
                    '-DOGS_CHECK_HEADER_COMPILATION=ON ' +
                    '-DOGS_BUILD_UTILS=ON ' +
                    '-DOGS_BUILD_GUI=ON ' +
                    '-DOGS_USE_PYTHON=ON ' +
                    '-DBUILD_SHARED_LIBS=ON '
                  dir = 'build-check-header'
                }
              }
              catch(err) {
                sh 'cat build-check-header/CMakeFiles/CMakeError.log'
                unstable('check-header failed!')
              }
            }
          }
        }
        // ************************* Tests-Large *******************************
        stage('Tests-Large') {
          when {
            beforeAgent true
            allOf {
              expression { return params.master_jobs && (stage_required.build || stage_required.full) }
              environment name: 'JOB_NAME', value: 'ufz/ogs/master'
            }
          }
          agent {
            dockerfile {
              filename 'Dockerfile.gcc.full'
              dir 'scripts/docker'
              label 'singularity1 || envinf1'
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
            }
            cleanup {
              dir('build') { deleteDir() }
            }
          }
        }
        // ************************** Sanitizer ********************************
        stage('Sanitizer') {
          when {
            beforeAgent true
            allOf {
              expression { return params.master_jobs && (stage_required.build || stage_required.full)}
              environment name: 'JOB_NAME', value: 'ufz/ogs/master'
            }
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
          environment {
            UBSAN_OPTIONS = 'print_stacktrace=1'
            LSAN_OPTIONS = "suppressions=$WORKSPACE/scripts/test/leak_sanitizer.suppressions"
          }
          steps {
            script {
              sh 'git submodule sync && git submodule update'
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
            allOf {
              expression { return params.master_jobs }
              environment name: 'JOB_NAME', value: 'ufz/ogs/master'
            }
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
                  --cmake ' -DOGS_USE_PYTHON=ON -DOGS_BUILD_UTILS=ON'
                python ThirdParty/container-maker/ogscm/cli.py -B -C -R \
                  -j $NUM_THREADS --ogs . --pm system --cvode \
                  --ompi 4.0.1
              '''.stripIndent()
            }
          }
          post {
            success {
              archiveArtifacts('_out/images/*.sif')
            }
            cleanup {
              dir('_out') { deleteDir() }
            }
          }
        }
        // ********************* Push Docker Images ***************************
        stage('Push Docker Images') {
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
              dir('scripts/docker') {
                def gccImage = docker.build("ogs6/gcc:latest", "-f Dockerfile.gcc.full .")
                def gccGuiImage = docker.build("ogs6/gcc:gui", "-f Dockerfile.gcc.gui .")
                def clangImage = docker.build("ogs6/clang:latest", "-f Dockerfile.clang.full .")
                withCredentials([usernamePassword(credentialsId: 'docker-hub-credentials',
                  passwordVariable: 'pw', usernameVariable: 'docker_user')]) {
                  sh 'echo $pw | docker login -u $docker_user --password-stdin'
                  retry(3) {
                    gccImage.push()
                    gccGuiImage.push()
                    clangImage.push()
                  }
                }
              }
            }
          }
        }
        // *************************** Web *************************************
        stage('Web') {
          agent {
            dockerfile {
              filename 'Dockerfile.web'
              dir 'scripts/docker'
              label 'docker'
              args '--entrypoint='
            }
          }
          when {
            beforeAgent true
            environment name: 'JOB_NAME', value: 'ufz/ogs/master'
          }
          steps {
            dir('web') {
              sh 'urlchecker check --retry-count 5 --file-types .pandoc,.md --white-listed-files releases/* --white-listed-urls https://jenkins.opengeosys.org/job/ufz/job/ogs-container-maker/job/master/build,http://yourproxy.example.com,https://apt.kitware.com/ubuntu/,https://github.com/YOUR-USERNAME/ogs,https://jenkins.opengeosys.org/github-webhook/,http://localhost:1313,https://github.com/ufz/ogs/pull/\\$1,http://www.opengeosys.org/images/xsd/OpenGeoSysXXX.xsd,https://\\`-protocol content'
            }
          }
        }
        // *************************** Post ************************************
        stage('Post') {
          agent { label "master"}
          when { buildingTag() }
          steps {
            script {
              currentBuild.keepLog(true)
              currentBuild.displayName = tag
            }
          }
        }
      } // end parallel
    } // end stage Build
  }
}
