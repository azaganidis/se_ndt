get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${SELF_DIR}/se_ndt.cmake)
#set(se_ndt_LIBRARIES @se_ndt_EXPORT_LIBRARIES@)
set(se_ndt_LIBRARIES @se_ndt_LIBRARIES@)
