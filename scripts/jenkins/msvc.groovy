defaultCMakeOptions = '-DCMAKE_BUILD_TYPE=Release -DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System -DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'

node('win && conan')
{
    stage 'Checkout (Win)'
    dir('ogs') { checkout scm }

    withEnv(getEnv()) {
        stage 'Configure (Win)'
        configure 'build', '', 'Ninja',
            '-u -s build_type=Release -s compiler="Visual Studio" -s compiler.version=12 -s arch=x86_64'

        stage 'CLI (Win)'
        build 'build', 'package'

        stage 'Test (Win)'
        build 'build', 'tests'

        stage 'Data Explorer (Win)'
        configure 'build', '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF',
            'Ninja', '-u -s build_type=Release -s compiler="Visual Studio" -s compiler.version=12 -s arch=x86_64',
            true
        build 'build', 'package'
    }

    if (env.BRANCH_NAME == 'master' || env.BRANCH_NAME.contains('release') ) {
        stage 'Release (Win)'
        archive 'build/*.zip'
    }

    stage 'Post (Win)'
    publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/msvc-log-parser.rules'
}

// *** Helper functions ***
def getEnv()
{
    if (env.NODE_NAME == 'visserv3')
        qtdir = 'C:\\libs\\qt\\4.8\\msvc2013-x64'
    if (env.NODE_NAME == 'win1')
        qtdir = 'C:\\libs\\qt-4.8.7-x64-msvc2013\\qt-4.8.7-x64-msvc2013'

    return [
        "QTDIR=${qtdir}",
        'Path=$Path;$QTDIR\\bin',
        'CONAN_CMAKE_GENERATOR=Ninja'
    ]
}


def configure(buildDir, cmakeOptions, generator, conan_args=null, keepBuildDir=false) {
    if (keepBuildDir == false)
        bat("""rd /S /Q ${buildDir}
               mkdir ${buildDir}""".stripIndent())
    if (conan_args != null)
        bat("""cd ${buildDir}
               conan install ../ogs ${conan_args}""".stripIndent())
    bat """set path=%path:\"=%
           call "%vs120comntools%..\\..\\VC\\vcvarsall.bat" x86_amd64
           cd ${buildDir}
           cmake ../ogs -G "${generator}" ${defaultCMakeOptions} ${cmakeOptions}"""
}

def build(buildDir, target=null) {
    targetString = ""
    if (target != null)
        targetString = "--target ${target}"
    bat("""set path=%path:\"=%
           call "%vs120comntools%..\\..\\VC\\vcvarsall.bat" x86_amd64
           cd ${buildDir}
           cmake --build . --config Release ${targetString}""".stripIndent())
}

def deploy(files) {
    step([$class: 'S3BucketPublisher', dontWaitForConcurrentBuildCompletion: true, entries: [[bucket: 'opengeosys', excludedFile: '', flatten: true, gzipFiles: false, managedArtifacts: true, noUploadOnFailure: true, selectedRegion: 'eu-central-1', sourceFile: "${files}", storageClass: 'STANDARD', uploadFromSlave: true, useServerSideEncryption: false]], profileName: 'S3 UFZ', userMetadata: []])
}

def publishTestReports(ctestPattern, gtestPattern, parseRulefile) {
    step([$class: 'XUnitPublisher', testTimeMargin: '3000', thresholdMode: 1,
        thresholds: [
            [$class: 'FailedThreshold', failureNewThreshold: '', failureThreshold: '', unstableNewThreshold: '', unstableThreshold: ''],
            [$class: 'SkippedThreshold', failureNewThreshold: '', failureThreshold: '', unstableNewThreshold: '', unstableThreshold: '']],
        tools: [
            [$class: 'GoogleTestType', deleteOutputFiles: true, failIfNotNew: true, pattern: "${gtestPattern}", skipNoTestFiles: false, stopProcessingIfError: true]]
    ])

    step([$class: 'LogParserPublisher', failBuildOnError: true, unstableOnWarning: false,
            projectRulePath: "${parseRulefile}", useProjectRule: true])
}
// *** End helper functions ***
