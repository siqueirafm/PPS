add_custom_target(doc
                  COMMAND doxygen "${CMAKE_CURRENT_SOURCE_DIR}/config.dox"
		  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		  VERBATIM
		  COMMENT "Generate documentation for the PPS project")

