
SET( display_sources
  cdcl_display.h                    cdcl_display.txx
)

# vgui and VTK for visualization (not required)
INCLUDE_DIRECTORIES( ${VXL_VGUI_INCLUDE_DIR} )

AUX_SOURCE_DIRECTORY(Templates display_sources)

# dependencies: 
# vgui and VTK for visualization (not required for cdcl)
INCLUDE_DIRECTORIES( ${VXL_VGUI_INCLUDE_DIR} )
IF( VXL_VGUI_FOUND )
  ADD_LIBRARY( cdcl_display ${display_sources} )
  IF( EXISTS ${VTK_USE_FILE} )
    TARGET_LINK_LIBRARIES( cdcl_display vgui vtkIO vtkCommon )
  ELSE( EXISTS ${VTK_USE_FILE} )
    TARGET_LINK_LIBRARIES( cdcl_display vgui )
  ENDIF( EXISTS ${VTK_USE_FILE} )
ENDIF( VXL_VGUI_FOUND )


