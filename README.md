# FilmDosimetryOptimization

This respository contains the core methodology for estimating EBT2/EBT3 doses
using multi-channel optimization. The metodology described at my [publication](https://www.sciencedirect.com/science/article/pii/S0010465513000805?via%3Dihub)
was ported to C++ using OpenCV.

More details at my doctorate's degree
[thesis](http://antigo.nuclear.ufrj.br/DScTeses/teses2014/Tese_Victor_Alves.pdf) (in Portuguese with english abstract)
That research aimed to develop robust optimization techniques on multi-channel film dosimetry, then to evaluate absorbed dose to water and its uncertainty.

The overall methodology is implemented on [Film2Dose software](https://conferences.iaea.org/indico/event/108/session/94/contribution/185.pdf).

## Getting Started

Check main.cxx file.

    Film2DoseOptimizaion.exe filePath.tiff

### Prerequisites

  * C++ compiler
  * C++11 support
  * OpenCV

## Installing
  * Tested on Windows 7/10 - 64 bits , Qt5.10.0, MSVC-2017 64 bits.
  * Lubuntu linux 16.04 + KDE plasma - 64 bits, Qt5.10.0, LLVM 5.0.1/Clang++

## Built With

* [CMAKE](https://cmake.org/)
* [QT-Creator](https://en.wikipedia.org/wiki/Qt_Creator)
* [Ninja-build](https://ninja-build.org/)

## Contributing

Any bug fixes or improvements are welcome.

## Versioning

[SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

## Author
    Victor Gabriel Leandro Alves, D.Sc.
    Copyright 2010-2018

## License

This project is distributed under the BSD-style licence

## Acknowledgments
* The best free C++ IDE if have found: [QT-creator Open Source](https://www.qt.io/download-qt-for-application-development)
* John Purcell's [free C++ course](https://www.udemy.com/free-learn-c-tutorial-beginners/) and [advanced C++ course](https://www.udemy.com/learn-advanced-c-programming/)
