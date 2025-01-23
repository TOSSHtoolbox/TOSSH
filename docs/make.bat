@ECHO OFF

pushd %~dp0

REM Command file for Sphinx documentation

if "%SPHINXBUILD%" == "" (
	set SPHINXBUILD=sphinx-build
)
set SOURCEDIR=.
set BUILDDIR=_build

if "%1" == "" goto help

%SPHINXBUILD% >NUL 2>NUL
if errorlevel 9009 (
	echo.
	echo.The 'sphinx-build' command was not found. Make sure you have Sphinx
	echo.installed, then set the SPHINXBUILD environment variable to point
	echo.to the full path of the 'sphinx-build' executable. Alternatively you
	echo.may add the Sphinx directory to PATH.
	echo.
	echo.If you don't have Sphinx installed, grab it from
	echo.https://www.sphinx-doc.org/
	exit /b 1
)

%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
goto end

:help
%SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%

:end
popd

:: Building the Documentation Locally before pushing to GitHub

:: (A) [Recommended] Steps to Enable Real-Time Updates

:: (1) Install sphinx, sphinx-rtd-theme, and sphinx-autobuild in your conda environment
:: (2) Navigate to the docs directory
:: (3) Run the make.bat file with the html option: make.bat html
:: (3)  Start a local HTTP Server ``` sphinx-autobuild . _build/html```. This starts a server at http://127.0.0.1:8000

:: (B) Build it static

:: (1) Install sphinx and phinx-rtd-theme in your conda environment
:: (2) Navigate to the docs directory
:: (3) Run the make.bat file with the html option: make.bat html
:: (4) Navigate to the _build/html directory
:: (5) Start a local HTTP Server ``` python -m http.server ```. This starts a server at http://localhost:8000
