node('docker')
{
	// Checks out into subdirectory ogs
	stage 'Checkout'
	checkout([$class: 'GitSCM',
		branches: [[name: '*/master']],
		doGenerateSubmoduleConfigurations: false,
		extensions:
			[[$class: 'RelativeTargetDirectory', relativeTargetDir: 'ogs']],
		submoduleCfg: [],
		userRemoteConfigs:
			[[url: 'https://github.com/ufz/ogs']]])


	// Multiple configurations are build in parallel
	parallel linux: {
		docker.image('ogs6/gcc-ogs-cli').inside
		{
			build 'build', '', 'package'

			stage 'Test'
			sh '''cd build
			      rm -rf tests/testrunner.xml
			      bin/testrunner --gtest_output=xml:./tests/testrunner.xml'''
		}
	},

	linux_gui: {
		docker.image('ogs6/gcc-ogs-gui').inside {
			build 'build_gui', '-DOGS_BUILD_GUI=ON -DOGS_BUILD_TESTS=OFF -DOGS_BUILD_CLI=OFF', 'package'
		}
	}

//	windows: {
//		docker.image('ogs6/mingw-ogs-gui').inside
//		{
//			build 'build_win', '-DCMAKE_TOOLCHAIN_FILE=$CMAKE_TOOLCHAIN_FILE', 'package'
//
//			stage 'Test'
//			sh '''cd build_win
//						rm -rf tests/testrunner.xml
//						wine bin/testrunner.exe --gtest_output=xml:./tests/testrunner.xml'''
//		}
//	},
//
//	windows_gui: {
//		docker.image('ogs6/mingw-ogs-gui').inside
//		{
//			build 'build_win_gui',
//						'-DCMAKE_TOOLCHAIN_FILE=$CMAKE_TOOLCHAIN_FILE -DOGS_BUILD_GUI=ON -DOGS_BUILD_TESTS=OFF -DOGS_BUILD_CLI=OFF',
//						'package'
//		}
//	}

	// end parallel

	step([$class: 'JUnitResultArchiver',
		testResults: 'build/tests/testrunner.xml,build_win/tests/testrunner.xml'])

	archive 'build*/*.tar.gz,build_win*/*.zip'
} // end node


def build(buildDir, cmakeOptions, target) {
	sh "rm -rf ${buildDir} && mkdir ${buildDir}"

	stage 'Configure'
	sh "cd ${buildDir} && cmake ../ogs ${cmakeOptions}"

	stage 'Build'
	sh "cd ${buildDir} && make -j 2 ${target}"
}
