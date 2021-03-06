/**
*
* BEMoDA v1.0 - BiowaivEr aid for Model Dependent-Independent Approach script for in-vitro dissolution profile comparison
Copyright (C) 2017 Jakub Szlęk, Aleksander Mendyk

Authors: 
Jakub Szlęk, Aleksander Mendyk

Affiliation: 
Jagiellonian University Medical College,
Faculty of Pharmacy,
Department of Pharmaceucial Technology and Biopharmaceutics,
Medyczna 9 st.,
30-688 Kraków
Poland

Bugs, issues, please e-mail to maintainer
Jakub Szlęk: j.szlek@uj.edu.pl

Copyright (C) 2017 Jakub Szlęk, Aleksander Mendyk

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General 
Public License as published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.
If not, see <http://www.gnu.org/licenses/>.

========================
REQUIREMENTS
========================
1) R environment version >= 3.3.0 (for download check https://cran.r-project.org/)
2) Additional R packages: dplyr, optimx, ggplot2, nloptr, GenSA, MASS

Note: Please ensure that nlopt libraries (https://nlopt.readthedocs.io/en/latest/) are installed on your system.
If you encounter troubles installing nlopt please read the StackOverflow question https://stackoverflow.com/questions/29716857/installing-nloptr-on-linux

========================
INSTALLATION AND RUNNING
========================

1) Download the tar.gz archive from the https://sourceforge.net/projects/bemoda/
2) Unpack downloaded package (on linux-like systems use tar xzvf BEMoDA.tar.gz -C /path/to/my/BEMoDA/folder)
3) Run the konsole/terminal or cmd.exe
4) Navigate the unpacked (on linux-like systems use cd /path/to/my/BEMoDA/folder)
5) On Windows run the run_R.bat file or on linux-like systems run run_R.sh (make sure the *.bat or *.sh files are executable)
