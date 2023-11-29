
pybind11_add_module(mcrl2py mcrl2py.cpp)

target_link_libraries(mcrl2py PRIVATE mcrl2_process mcrl2_lps)

add_test(NAME mcrl2py_test
  COMMAND
  ${CMAKE_COMMAND} -E env
    PYTHONPATH=$<TARGET_FILE_DIR:mcrl2py>
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/mcrl2py_test.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
