def post = new ogs.post()

stage('Docset') {
    sh("""
    rm -rf build; mkdir build; cd build
    cmake -DDOCS_GENERATE_DOCSET=ON ../ogs
    cmake --build . --config Release --target doc

    cd docs
    tar --exclude='.DS_Store' -cvzf ogs6.tgz ogs6.docset
    """.stripIndent())
    archiveArtifacts 'build/docs/ogs6.tgz, build/docs/ogs6.xml'

    post.cleanup()
}
