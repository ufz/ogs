if(WIN32)
    set(_envrc_file_ending ".ps1")
    set(_envrc_content
        "$Env:PYTHONPATH += \"\;${PROJECT_BINARY_DIR}/site-packages\""
        "$Env:PATH += \"\;${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}\""
        "$Env:OGS_USE_PATH = \"1\""
    )
    if(OGS_USE_PIP)
        list(APPEND _envrc_content
            "$Env:UV_PROJECT = \"${PROJECT_SOURCE_DIR}/Tests/Data\""
            "$Env:UV_PROJECT_ENVIRONMENT = \"${PROJECT_BINARY_DIR}/.venv\""
            "$Env:UV_FROZEN = \"1\""
        )
    endif()
    if(OGS_USE_MKL)
        list(APPEND _envrc_content
            "$Env:PATH += \"\;C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin\""
            "$Env:PATH += \"\;C:/Program Files (x86)/Intel/oneAPI/mkl/latest/bin\""
            "Invoke-BatchFile \"C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\latest\\env\\vars.bat\""
            "Invoke-BatchFile \"C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\latest\\env\\vars.bat\""
        )
    endif()
else()
    set(_envrc_content
        "[ -d \"${PROJECT_BINARY_DIR}/.venv\" ] && source \"${PROJECT_BINARY_DIR}/.venv/bin/activate\""
        "export PATH=\"${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}:$PATH\""
        "export OGS_USE_PATH=1"
        "export PYTHONPATH=\"${PROJECT_BINARY_DIR}/site-packages:$PYTHONPATH\""
    )
    if(OGS_USE_PIP)
        list(APPEND _envrc_content
            "export UV_PROJECT=\"${PROJECT_SOURCE_DIR}/Tests/Data\""
            "export UV_PROJECT_ENVIRONMENT=\"${PROJECT_BINARY_DIR}/.venv\""
            "export UV_FROZEN=1"
        )
    endif()
    if(TFEL_WITH_PYTHON)
        set(_envrc_content "${_envrc_content}"
                           "export PYTHONPATH=${TFEL_WITH_PYTHON}:$PYTHONPATH"
        )
    endif()
endif()
string(JOIN "\n" _envrc_content_text ${_envrc_content})
file(CONFIGURE OUTPUT .envrc${_envrc_file_ending} CONTENT "${_envrc_content_text}")
