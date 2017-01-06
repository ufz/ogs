macro(InstallXmlSchemaFiles GLOB_EXPRESSION)
    file(GLOB XSD_FILES . ${GLOB_EXPRESSION})
    install(FILES ${XSD_FILES} DESTINATION bin COMPONENT ogs_cli)
endmacro()
