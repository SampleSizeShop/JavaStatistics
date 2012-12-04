%%
%% The `glimmpse-report.cls' is used for formatting 
%% GLIMMPSE validation reports.
%%
%% Copyright (C) 2010 Regents of the University of Colorado.  
%%
%% This program is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public License
%% as published by the Free Software Foundation; either version 2
%% of the License, or (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
%%
%% Author: Sarah Kreidler
%%
\def\fileversion{1.0}
\def\filename{glimmpse-report}
\def\filedate{2012/09/06}
%%
%% Package 'glimmpse-report' to use with LaTeX2e 
%% License: GPL-2
%% Copyright: (C) University of Colorado Denver
%% Author: Sarah Kreidler
%% Please report errors to sarah.kreidler@ucdenver.edu
%%
%% Based on the package 'jss' of the Journal of Statistical Software
%% See http://www.jstatsoft.org/ for details
%%
%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{glimmpse-report}[\filedate\space\fileversion\space Glimmpse Validation Report Class by Sarah Kreidler]
\ProcessOptions

%% load article defaults
\LoadClass[12pt,english]{article}

%% required packages
\RequirePackage{savesym,amsmath,graphicx,color,hyperref,geometry,fancyhdr,extramarks,lastpage,chngpage}
\savesymbol{iint}
\usepackage[T1]{fontenc}
\usepackage[latin9]{inputenc}
\usepackage{esint}
\usepackage[english]{babel}
%% long tables that span multiple pages
\renewcommand{\familydefault}{\sfdefault}
\usepackage{longtable,tabu}
\PassOptionsToPackage{normalem}{ulem}
\usepackage{ulem}
\usepackage{tabularx}

%% set margins
\geometry{verbose,tmargin=1in,bmargin=1in,lmargin=0.5in,rmargin=0.5in}
\setlength{\parskip}{\medskipamount}
\setlength{\parindent}{0pt}

%% definitions from calling document
\newcommand{\JavaStatisticsVersion}[1]{\def\@JavaStatisticsVersion{#1}}
\newcommand{\doctitle}[1]{\def\@doctitle{#1}}
\newcommand{\docauthor}[1]{\def\@docauthor{#1}}
\newcommand{\docdate}[1]{\def\@docdate{#1}}

%% header and footer setup
\pagestyle{fancy}

%% header
\lhead{\includegraphics[height=13pt]{glimmpse.png} GLIMMPSE Validation Report}
\chead{} 
\rhead{JavaStatistics Library Version \@JavaStatisticsVersion} %

%% footer
\lfoot{\includegraphics[height=13pt]{UCD_CSPH.png} \includegraphics[height=13pt]{UCD_CSPH_text.png}}                                                      
\cfoot{\includegraphics[height=13pt]{UFL_HOP.png}}                                                              
\rfoot{Page\ \thepage\ of\ \protect\pageref{LastPage}}                  %

%% header / footer lines
\renewcommand\headrulewidth{0.4pt}                                      %
\renewcommand\footrulewidth{0.4pt}                                      %

%% colors
\definecolor{Red}{rgb}{0.5,0,0}
\definecolor{Blue}{rgb}{0,0,0.5}

%% hyperlink formatting
\hypersetup{%
    hyperindex = {true},
    colorlinks = {true},
    linktocpage = {true},
    plainpages = {false},
    linkcolor = {Blue},
    citecolor = {Blue},
    urlcolor = {Blue},
    pdfstartview = {Fit},
    pdfpagemode = {UseOutlines},
    pdfview = {XYZ null null null}
}

\makeatletter

%% Add title information for 
%% validation reports
\AtBeginDocument{
  \setcounter{MaxMatrixCols}{100}
  
  {\huge \bfseries GLIMMPSE Validation Report: }
  
  {\LARGE \@doctitle}
  
  % Author
  \emph{Authors:} \@docauthor
  
  % Date
  \emph{Run Date:} \@docdate

}

%% bibliography
\usepackage[authoryear,round,longnamesfirst]{natbib}
\bibpunct{(}{)}{;}{a}{}{,}
\bibliographystyle{glimmpse}
%% paragraphs
\setlength{\parskip}{0.7ex plus0.1ex minus0.1ex}
\setlength{\parindent}{0em}

%% sections, subsections, and subsubsections
\newlength{\preXLskip}
\newlength{\preLskip}
\newlength{\preMskip}
\newlength{\preSskip}
\newlength{\postMskip}
\newlength{\postSskip}
\setlength{\preXLskip}{1.8\baselineskip plus 0.5ex minus 0ex}
\setlength{\preLskip}{1.5\baselineskip plus 0.3ex minus 0ex}
\setlength{\preMskip}{1\baselineskip plus 0.2ex minus 0ex}
\setlength{\preSskip}{.8\baselineskip plus 0.2ex minus 0ex}
\setlength{\postMskip}{.5\baselineskip plus 0ex minus 0.1ex}
\setlength{\postSskip}{.3\baselineskip plus 0ex minus 0.1ex}

\newcommand{\jsssec}[2][default]{
  \vskip \preXLskip%
  \pdfbookmark[1]{#1}{Section.\thesection.#1}%
  \refstepcounter{section}%
  \centerline{\textbf{\Large \thesection. #2}} \nopagebreak
  \vskip \postMskip \nopagebreak
  }
\newcommand{\jsssecnn}[1]{\vskip \preXLskip%
  \centerline{\textbf{\Large #1}} \nopagebreak
  \vskip \postMskip \nopagebreak}

\newcommand{\jsssubsec}[2][default]{\vskip \preMskip%
  \pdfbookmark[2]{#1}{Subsection.\thesubsection.#1}%
  \refstepcounter{subsection}%
  \textbf{\large \thesubsection. #2} \nopagebreak
  \vskip \postSskip \nopagebreak}
\newcommand{\jsssubsecnn}[1]{\vskip \preMskip%
  \textbf{\large #1} \nopagebreak
  \vskip \postSskip \nopagebreak}

\newcommand{\jsssubsubsec}[2][default]{\vskip \preSskip%
  \pdfbookmark[3]{#1}{Subsubsection.\thesubsubsection.#1}%
  \refstepcounter{subsubsection}%
  {\large \textit{\thesubsubsection. #2}} \nopagebreak
  \vskip \postSskip \nopagebreak}
\newcommand{\jsssubsubsecnn}[1]{\vskip \preSskip%
  {\textit{\large #1}} \nopagebreak
  \vskip \postSskip \nopagebreak}

\newcommand{\jsssimplesec}[2][default]{\vskip \preLskip%
  \pdfbookmark[1]{#1}{Section.\thesection.#1}%
  \refstepcounter{section}%
  \textbf{\large #1} \nopagebreak
  \vskip \postSskip \nopagebreak}
\newcommand{\jsssimplesecnn}[1]{\vskip \preLskip%
  \textbf{\large #1} \nopagebreak
  \vskip \postSskip \nopagebreak}

\renewcommand{\section}{\secdef \jsssec \jsssecnn}
\renewcommand{\subsection}{\secdef \jsssubsec \jsssubsecnn}
\renewcommand{\subsubsection}{\secdef \jsssubsubsec \jsssubsubsecnn}

%%
%% End of file `glimmpse.cls'.
