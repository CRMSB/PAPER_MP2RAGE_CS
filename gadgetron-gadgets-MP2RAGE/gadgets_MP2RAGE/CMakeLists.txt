set( gadgetron_rmsb_header_files ExtractMP2RAGEGadget.h

  	)

set( gadgetron_rmsb_src_files   ExtractMP2RAGEGadget.cpp

        )

set( gadgetron_rmsb_config_files
	config/MP2RAGE_CS_Bucket.xml
	config/MP2RAGE_CS_Bucket_DUMP.xml
)

add_library(gadgetron_rmsb SHARED  
${gadgetron_rmsb_header_files} 
${gadgetron_rmsb_src_files}
${gadgetron_rmsb_config_files}
)
target_link_libraries(gadgetron_rmsb Gadgetron::gadgetron_core Gadgetron::gadgetron_toolbox_cpucore_math)


install(TARGETS gadgetron_rmsb
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES ${gadgetron_rmsb_config_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)