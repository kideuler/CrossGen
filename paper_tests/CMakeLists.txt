# The experiments of Sec. 6 of the paper, E1 to E5.
#
# One library of shared machinery (domains, the three methods behind one
# interface, the Sec. 5 metrics, the CSV/console reporting) and one executable
# per experiment. Each executable writes its tables to --out (default
# ./results) and exits non-zero if a check it makes fails, so they can be run
# by hand for the numbers or by ctest as regression tests.

# The Qt viewer turns CMAKE_AUTOMOC on globally and it is inherited here, where
# nothing has a Q_OBJECT in it. Off, so the build does not run moc over every
# file for nothing.
set(CMAKE_AUTOMOC OFF)

add_library(PaperCommon STATIC
    common/Domains.hxx
    common/Domains.cxx
    common/Metrics.hxx
    common/Metrics.cxx
    common/Methods.hxx
    common/Methods.cxx
    common/Report.hxx
    common/Corpus.hxx
    common/Layout.hxx
    common/Layout.cxx
)
target_include_directories(PaperCommon PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_link_libraries(PaperCommon PUBLIC PolyVector)

foreach(EXP E1_Verification E2_Canonical E3_Benchmark E4_Layout E5_MultiMaterial MakeDomains)
    add_executable(Paper${EXP} ${EXP}.cxx)
    target_link_libraries(Paper${EXP} PRIVATE PaperCommon)
    set_target_properties(Paper${EXP} PROPERTIES OUTPUT_NAME ${EXP})
endforeach()

# Where the corpus lives, so an experiment does not have to be told.
foreach(EXP E1_Verification E2_Canonical E3_Benchmark E4_Layout E5_MultiMaterial MakeDomains)
    target_compile_definitions(Paper${EXP} PRIVATE
        PAPER_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../data/meshes")
endforeach()

# E1 and E2 are quick and are run by ctest; E3 to E5 sweep the whole corpus and
# are not, because a full sweep is minutes to hours and belongs in a run of its
# own. `make paper` is that run.
add_test(NAME Paper_E1_Verification COMMAND E1_Verification --out ${CMAKE_CURRENT_BINARY_DIR}/results --levels 3)
set_tests_properties(Paper_E1_Verification PROPERTIES TIMEOUT 3600)
add_test(NAME Paper_E2_Canonical COMMAND E2_Canonical --out ${CMAKE_CURRENT_BINARY_DIR}/results --quick)
set_tests_properties(Paper_E2_Canonical PROPERTIES TIMEOUT 3600)

add_custom_target(paper
    COMMAND $<TARGET_FILE:PaperE1_Verification> --out ${CMAKE_CURRENT_BINARY_DIR}/results
    COMMAND $<TARGET_FILE:PaperE2_Canonical>    --out ${CMAKE_CURRENT_BINARY_DIR}/results
    COMMAND $<TARGET_FILE:PaperE3_Benchmark>    --out ${CMAKE_CURRENT_BINARY_DIR}/results
    COMMAND $<TARGET_FILE:PaperE4_Layout>       --out ${CMAKE_CURRENT_BINARY_DIR}/results
    COMMAND $<TARGET_FILE:PaperE5_MultiMaterial> --out ${CMAKE_CURRENT_BINARY_DIR}/results
    DEPENDS PaperE1_Verification PaperE2_Canonical PaperE3_Benchmark PaperE4_Layout PaperE5_MultiMaterial
    COMMENT "Running the paper's Tier-1 experiments E1-E5"
    VERBATIM
)
