function(check_path_and_inform MODE LIB PATH NAME)
    if(NOT EXISTS ${PATH})
        message(${MODE} "${PATH} does not exist, specify the path to a working\
    ${LIB} installation by -D${NAME}=...")
    else()
        message(STATUS "Found ${LIB}: ${PATH}")
    endif()
endfunction()