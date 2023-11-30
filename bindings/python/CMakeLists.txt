
pybind11_add_module(mcrl2py mcrl2_bindings.cpp)


target_link_libraries(mcrl2py PRIVATE mcrl2_process mcrl2_lps)