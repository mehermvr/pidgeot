include(FetchContent)
FetchContent_Declare(
    cli11
    QUIET
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG f4d0731
)
FetchContent_MakeAvailable(cli11)
