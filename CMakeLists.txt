Project(GBPLSegmentation)

include_directories(./include)

add_executable(tttGBPLSamplerTest src/tttGBPLSamplerTest.cpp)

target_link_libraries(tttGBPLSamplerTest ${ITK_LIBRARIES})

add_executable(tttGBPLInferenceTest src/tttGBPLInferenceTest.cpp)

target_link_libraries(tttGBPLInferenceTest ${ITK_LIBRARIES})
