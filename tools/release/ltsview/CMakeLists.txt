mcrl2_add_gui_tool(ltsview
  MENUNAME "LTSView"
  DESCRIPTION "Visualisation tool for large state spaces (mCRL2 toolset)"
  ICON "mcrl2-green"
  SOURCES
    ./source/visualizer/
    ./source/visualizer/generic_visualizer/
    ./source/visualizer/generic_visualizer/cameras/
    ./source/visualizer/ltsview_visualizer/
    ./source/
    ./source/ui/
    # Qt wants .ui files in the sources variable; include by directory.
    ./include/ui/ui_files/
  DEPENDS
    mcrl2_lts
    mcrl2_gui
    Qt6::Core
    Qt6::Gui
    Qt6::OpenGL
    Qt6::OpenGLWidgets
    Qt6::Widgets
    Qt6::Xml
    OpenGL::GL
    OpenGL::GLU
    tr
)
