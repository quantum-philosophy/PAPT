cmake_minimum_required (VERSION 2.6)

project(rota)

set(VENDOR_DIR ${PROJECT_SOURCE_DIR}/vendor)

set(HOTSPOT_VERSION 5.02)

add_subdirectory(msrc)
add_subdirectory(vendor)

add_custom_target(tags
	COMMAND cd ${PROJECT_BINARY_DIR} && ctags -R ../*
)
