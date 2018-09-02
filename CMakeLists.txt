cmake_minimum_required(VERSION 3.10)
project(program C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(program
        d_util.c
        d_util.h
        genpop.c
        genpop.h
        main.c
        model_Autofecondation.c
        model_CycleDemographique.c
        model_SexRatio.c
        model_Variance.c
        model_WrightFisher.c)
