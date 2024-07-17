# SWIP
Surface-Wave dispersion Inversion and Profiling software
Installation instructions and manual are stored in SWIP/docs/manual

# What's new
- Works under Windows 10 or higher with the Windows Subsystem for Linux
- Includes a first arrival picking tool for Seismic Refraction Tomography (saves in .sgt format for [pyGIMLi](https://www.pygimli.org/))

# Installation

1. Install MATLAB
2. Download SWIP at https://github.com/SWIPdev/SWIP/releases and add the content of the archive
to your MATLAB path
3. (Windows user only) Install Windows Subsystem for Linux (https://pureinfotech.com/install-wsl-windows-11/)
4. (Mac OS user only) Install Xcode, Homebrew and gnu-tar
5. Install Seismic Unix (https://github.com/JohnWStockwellJr/SeisUnix)
6. Install Geospy (http://www.geopsy.org/download.php)
7. Install SU extra binaries located in SWIP/src (for Mac users, install seg2segy located in Seismic Unix's third party)
8. (Optional) Install ImageMagick, PDFjam, and pdfcrop
