
function( add_example module name )
	add_executable( ${name} ${name}.cpp )
	target_link_libraries( ${name} ${module} )
	set_property( TARGET ${name}
		PROPERTY RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}" )
	set_property( TARGET ${name} PROPERTY FOLDER examples )
endfunction()

add_example(gpu_optimize_storm_drift test_optimize_storm_drift)
