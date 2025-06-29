@echo off
setlocal enabledelayedexpansion

:: Navigate five directories up to locate version.json
set "VERSION_FILE=%~dp0..\..\..\..\..\version.json"

:: Debug: Check if the file exists
if not exist "%VERSION_FILE%" (
    echo ERROR: version.json not found at %VERSION_FILE%
    exit /b 1
)

:: Read the version from version.json
set "VERSION=Unknown"
for /f "usebackq tokens=2 delims=:," %%A in ("%VERSION_FILE%") do (
    set "VERSION=%%A"
    set "VERSION=!VERSION:~2,-1!"
)

:: Debug: Check if a valid version was extracted
if "%VERSION%"=="Unknown" (
    echo ERROR: Failed to extract version from version.json
    exit /b 1
)

:: Define the path for version_define.h (same directory as the script)
set "OUTPUT_FILE=%~dp0version_define.h"

:: Create version_define.h
(
    echo #ifndef VERSION_DEFINE_H
    echo #define VERSION_DEFINE_H
    echo #define PROJECT_VERSION "!VERSION!"
    echo #endif
) > "%OUTPUT_FILE%"

:: Debug: Confirm file creation
if not exist "%OUTPUT_FILE%" (
    echo ERROR: Failed to create version_define.h
    exit /b 1
)

echo Created version_define.h with PROJECT_VERSION=!VERSION!
exit /b 0
