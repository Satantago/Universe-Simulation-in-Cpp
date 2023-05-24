#perf and behavior tests
add_executable(performances performances.cxx)
add_executable(test_pt test_point.cxx)
add_executable(test_uni test_univers.cxx)
add_executable(test_calcul test_calcul.cxx)
add_executable(conditions_limites conditions_limites.cxx)

target_link_libraries(performances particules)
target_link_libraries(test_pt particules)
target_link_libraries(test_uni particules)
target_link_libraries(test_calcul particules)
target_link_libraries(conditions_limites particules)


#unitary tests
add_executable(vecteur_test vecteur_test.cxx)
add_executable(grid_test grid_test.cxx)
add_executable(particule_test particule_test.cxx)
add_executable(univers_test univers_test.cxx)

target_link_libraries(vecteur_test gtest_main)
target_link_libraries(vecteur_test particules)

target_link_libraries(grid_test gtest_main)
target_link_libraries(grid_test particules)

target_link_libraries(particule_test gtest_main)
target_link_libraries(particule_test particules)

target_link_libraries(univers_test gtest_main)
target_link_libraries(univers_test particules)

include(GoogleTest)
gtest_discover_tests(vecteur_test)
gtest_discover_tests(grid_test)
gtest_discover_tests(particule_test)
gtest_discover_tests(univers_test)
