defaultDockerArgs = '-v /home/core/.ccache:/usr/src/.ccache'
defaultCMakeOptions = '-DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System'

node('docker')
{
	step([$class: 'GitHubSetCommitStatusBuilder', statusMessage: [content: 'Started Jenkins gcc build']])

	stage 'Checkout'
	dir('ogs') { checkout scm }

	stage 'Build'
	docker.image('ogs6/gcc-base:latest')
		.inside(defaultDockerArgs) {
			build 'build', '', 'package tests ctest'
		}

	publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
		'ogs/scripts/jenkins/clang-log-parser.rules'

	if (env.BRANCH_NAME == 'master')
		archive 'build*/*.tar.gz'
}


def build(buildDir, cmakeOptions, target) {
	sh "rm -rf ${buildDir} && mkdir ${buildDir}"

	stage 'Configure'
	sh "cd ${buildDir} && cmake ../ogs ${defaultCMakeOptions} ${cmakeOptions}"

	stage 'Build'
	sh "cd ${buildDir} && make -j \$(nproc) ${target}"
}

def publishTestReports(ctestPattern, gtestPattern, parseRulefile) {
	step([$class: 'XUnitPublisher', testTimeMargin: '3000', thresholdMode: 1,
		thresholds: [
			[$class: 'FailedThreshold', failureNewThreshold: '', failureThreshold: '', unstableNewThreshold: '', unstableThreshold: ''],
			[$class: 'SkippedThreshold', failureNewThreshold: '', failureThreshold: '', unstableNewThreshold: '', unstableThreshold: '']],
		tools: [
			[$class: 'CTestType', deleteOutputFiles: true, failIfNotNew: true, pattern: "${ctestPattern}", skipNoTestFiles: false, stopProcessingIfError: true],
			[$class: 'GoogleTestType', deleteOutputFiles: true, failIfNotNew: true, pattern: "${gtestPattern}", skipNoTestFiles: false, stopProcessingIfError: true]]
	])

	step([$class: 'LogParserPublisher', failBuildOnError: true, unstableOnWarning: false,
			projectRulePath: "${parseRulefile}", useProjectRule: true])

	step([$class: 'GitHubCommitNotifier', resultOnFailure: 'FAILURE', statusMessage: [content: 'Finished Jenkins gcc build']])

}
