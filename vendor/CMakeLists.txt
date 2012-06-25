set(HOTSPOT_VERSION 5.02)

set(HOTSPOT_URL http://lava.cs.virginia.edu/HotSpot/grab/HotSpot-${HOTSPOT_VERSION}.tar.gz)
set(HOTSPOT_SRC ${CMAKE_CURRENT_SOURCE_DIR}/HotSpot-${HOTSPOT_VERSION}.tar.gz)
set(HOTSPOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/HotSpot-${HOTSPOT_VERSION})

add_custom_target(hotspot
	COMMAND ${CMAKE_BUILD_TOOL} -C ${HOTSPOT_DIR}
	DEPENDS ${HOTSPOT_DIR}
)

add_custom_command(
	TARGET hotspot
	POST_BUILD
	COMMAND ln -fs ${HOTSPOT_DIR}/hotspot ${CMAKE_CURRENT_BINARY_DIR}/hotspot
	COMMAND ln -sfn ${HOTSPOT_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/hotspot
)

add_custom_command(OUTPUT ${HOTSPOT_DIR}
	COMMAND tar -xzf ${HOTSPOT_SRC} -C ${VENDOR_DIR}
	DEPENDS ${HOTSPOT_SRC}
)

add_custom_command(OUTPUT ${HOTSPOT_SRC}
	COMMAND wget -O ${HOTSPOT_SRC} ${HOTSPOT_URL}
)
