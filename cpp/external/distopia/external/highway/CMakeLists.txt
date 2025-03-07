cmake_minimum_required(VERSION 3.14)
project(Highway)

include(FetchContent)

set(HWY_GIT_REPO "https://github.com/google/highway.git")
set(HWY_GIT_TAG "5c688f8e9a29984ef59d07e5b5dfd8c5b1d8d8c5")

FetchContent_Declare(
    hwy
    GIT_REPOSITORY ${HWY_GIT_REPO}
    GIT_TAG ${HWY_GIT_TAG}
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(hwy)
